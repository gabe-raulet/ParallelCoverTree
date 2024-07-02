template <class PointTraits, class Index>
TreeIndex<PointTraits, Index>::TreeIndex()
    : points(Comm::comm_self()) {}

template <class PointTraits, class Index>
TreeIndex<PointTraits, Index>::TreeIndex(const Comm& comm)
    : points(comm) {}

template <class PointTraits, class Index>
TreeIndex<PointTraits, Index>::TreeIndex(const TreeIndex& rhs)
    : base(rhs.base),
      points(rhs.points) {}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::swap(TreeIndex& rhs)
{
    ::swap(base, rhs.base);
    points.swap(rhs.points);
}

template <class PointTraits, class Index>
TreeIndex<PointTraits, Index>&
TreeIndex<PointTraits, Index>::operator=(const TreeIndex& rhs)
{
    TreeIndex tmp(rhs);
    tmp.swap(*this);
    return *this;
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::build(const PointVector& mypoints, Real cutoff, Real base)
{
    this->base = base;
    points.assign(mypoints);
    init(cutoff);
}

template <class PointTraits, class Index>
Index TreeIndex<PointTraits, Index>::add_vertex(Index point, Index parent)
{
    Index vtx_level;
    Index vtx = pt.size();

    pt.push_back(point);
    children.emplace_back();

    if (parent >= 0)
    {
        vtx_level = level[parent] + 1;
        children[parent].push_back(vtx);
    }
    else vtx_level = 0;

    nlevels = max(vtx_level+1, nlevels);
    level.push_back(vtx_level);

    return vtx;
}

template <class PointTraits, class Index>
Index TreeIndex<PointTraits, Index>::batch_new_vertex(Index point, Index parent)
{
    my_new_vtx_pt_ids.push_back(point);
    my_new_vtx_hub_ids.push_back(parent);

    return pt.size() + my_new_vtx_pt_ids.size();
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::add_batched_vertices()
{
    auto comm = points.getcomm();

    assert(my_new_vtx_pt_ids.size() == my_new_vtx_hub_ids.size());

    using Iter = typename IndexVector::iterator;

    Iter ptitr, ptend, hubitr;
    IndexVector new_vtx_pt_ids, new_vtx_hub_ids;

    if (!comm.is_distributed())
    {
        ptitr = my_new_vtx_pt_ids.begin();
        ptend = my_new_vtx_pt_ids.end();
        hubitr = my_new_vtx_hub_ids.begin();
    }
    else
    {
        comm.allgatherv(my_new_vtx_pt_ids, new_vtx_pt_ids);
        comm.allgatherv(my_new_vtx_hub_ids, new_vtx_hub_ids);

        ptitr = new_vtx_pt_ids.begin();
        ptend = new_vtx_pt_ids.end();
        hubitr = new_vtx_hub_ids.begin();
    }

    while (ptitr != ptend)
        add_vertex(*ptitr++, *hubitr++);

    my_new_vtx_pt_ids.clear();
    my_new_vtx_hub_ids.clear();
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::initialize_root_hub()
{
    auto comm = getcomm();
    Index mysize = points.getmysize();

    my_dists.resize(mysize);
    my_hub_vtx_ids.resize(mysize);
    my_hub_pt_ids.resize(mysize);

    Index root = add_vertex(0, -1);
    hub_chains.insert({root, {pt[root]}});

    auto distance = PointTraits::distance();
    Point root_pt = points[0];

    if (comm.is_distributed()) comm.bcast(root_pt, 0);

    max_radius = -1;

    for (Index i = 0; i < mysize; ++i)
    {
        my_dists[i] = distance(root_pt, points[i]);
        my_hub_vtx_ids[i] = root;
        my_hub_pt_ids[i] = pt[root];
        max_radius = max(max_radius, my_dists[i]);
    }

    if (comm.is_distributed()) comm.allreduce(max_radius, MPI_MAX);
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::compute_farthest_hub_pts()
{
    auto comm = points.getcomm();
    Index mysize = points.getmysize();
    Index myoffset = points.getmyoffset();

    DistPairMap my_argmaxes;
    transform(hub_chains.begin(), hub_chains.end(), inserter(my_argmaxes, my_argmaxes.end()),
            [](auto pair) { return make_pair(pair.first, make_pair(-1, -1.0)); });

    for (Index i = 0; i < mysize; ++i)
    {
        Index hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            auto& [id, dist] = my_argmaxes.find(hub_id)->second;

            if (my_dists[i] > dist)
            {
                id = i;
                dist = my_dists[i];
            }
        }
    }

    if (!comm.is_distributed())
    {
        farthest_hub_pts = ::move(my_argmaxes);
    }
    else
    {
        using ArgmaxPair = MPIEnv::ArgmaxPair<Real, Index>;
        using ArgmaxPairVector = vector<ArgmaxPair>;

        IndexVector hub_ids;
        ArgmaxPairVector my_argmax_pairs;

        hub_ids.reserve(my_argmaxes.size());
        my_argmax_pairs.reserve(my_argmaxes.size());

        for (const auto& [hub_id, distpair] : my_argmaxes)
        {
            const auto& [id, dist] = distpair;

            hub_ids.push_back(hub_id);
            my_argmax_pairs.emplace_back(id + myoffset, dist);
        }

        MPI_Op MPI_ARGMAX;
        MPI_Datatype MPI_ARGMAX_PAIR;

        ArgmaxPair::create_mpi_handlers(MPI_ARGMAX_PAIR, MPI_ARGMAX);

        MPI_Allreduce(MPI_IN_PLACE, my_argmax_pairs.data(), static_cast<int>(my_argmax_pairs.size()), MPI_ARGMAX_PAIR, MPI_ARGMAX, comm.getcomm());

        MPI_Op_free(&MPI_ARGMAX);
        MPI_Type_free(&MPI_ARGMAX_PAIR);

        farthest_hub_pts.clear();

        for (Index i = 0; const auto& [id, dist] : my_argmax_pairs)
        {
            farthest_hub_pts.insert({hub_ids[i++], {id, dist}});
        }
    }
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::update_hub_chains()
{
    split_chains.clear(); leaf_chains.clear();

    for (const auto& [hub_id, farthest_pt] : farthest_hub_pts)
    {
        Real cover = level_radius(level[hub_id]) / base;
        auto [farthest_pt_id, farthest_dist] = farthest_pt;
        farthest_dist /= max_radius;

        if (farthest_dist == 0)
        {
            hub_chains.erase(hub_id);
            leaf_chains.insert(hub_id);
        }
        else if (farthest_dist <= cover)
        {
            split_chains.push_back(hub_id);
        }
        else
        {
            hub_chains.find(hub_id)->second.push_back(farthest_pt_id);
        }
    }
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::process_leaf_chains()
{
    Index mysize = points.getmysize();
    Index myoffset = points.getmyoffset();

    if (!leaf_chains.empty())
    {
        for (Index i = 0; i < mysize; ++i)
        {
            Index hub_id = my_hub_vtx_ids[i];

            if (leaf_chains.find(hub_id) != leaf_chains.end())
            {
                batch_new_vertex(i + myoffset, hub_id);
                my_hub_vtx_ids[i] = my_hub_pt_ids[i] = -1;
                my_dists[i] = 0;
            }
        }
    }

    add_batched_vertices();
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::process_split_chains()
{
    using IndexMap = unordered_map<Index, Index>;

    Index mysize = points.getmysize();

    if (!split_chains.empty())
    {
        IndexMap hub_pt_id_updates;
        IndexVector new_hub_pts, new_hub_ids;
        vector<IndexVector*> chain_ptrs(split_chains.size(), nullptr);
        transform(split_chains.begin(), split_chains.end(), chain_ptrs.begin(),
                [&](Index hub_id) { return &hub_chains.find(hub_id)->second; });

        Index num_new_hub_pts = 0;
        for (auto& chain_ptr : chain_ptrs)
            num_new_hub_pts += chain_ptr->size();

        new_hub_pts.reserve(num_new_hub_pts);
        new_hub_ids.reserve(num_new_hub_pts);

        for (Index i = 0; i < chain_ptrs.size(); ++i)
        {
            new_hub_pts.insert(new_hub_pts.end(), chain_ptrs[i]->cbegin(), chain_ptrs[i]->cend());
            new_hub_ids.insert(new_hub_ids.end(), chain_ptrs[i]->size(), split_chains[i]);
            hub_chains.erase(split_chains[i]);
        }

        for (Index i = 0; i < new_hub_pts.size(); ++i)
        {
            Index vtx_id = add_vertex(new_hub_pts[i], new_hub_ids[i]);
            hub_chains.insert({vtx_id, {new_hub_pts[i]}});
            hub_pt_id_updates.insert({new_hub_pts[i], vtx_id});
        }

        for (Index i = 0; i < mysize; ++i)
        {
            Index closest_pt_id = my_hub_pt_ids[i];
            auto it = hub_pt_id_updates.find(closest_pt_id);

            if (it != hub_pt_id_updates.end())
                my_hub_vtx_ids[i] = it->second;
        }
    }
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::update_dists_and_pointers()
{
    auto distance = PointTraits::distance();
    auto comm = points.getcomm();
    Index mysize = points.getmysize();
    bool distributed = comm.is_distributed();
    PointMap last_chain_pt_map;

    if (distributed)
    {
        IndexVector my_last_chain_pt_ids;
        transform(hub_chains.begin(), hub_chains.end(), back_inserter(my_last_chain_pt_ids),
                  [](const auto& pair) { return pair.second.back(); });

        PointVector last_chain_pts;
        IndexVector last_chain_pt_ids;

        points.allgather(my_last_chain_pt_ids, last_chain_pts, last_chain_pt_ids);

        for (Index i = 0; i < last_chain_pts.size(); ++i)
        {
            last_chain_pt_map.insert({last_chain_pt_ids[i], last_chain_pts[i]});
        }
    }

    for (Index i = 0; i < mysize; ++i)
    {
        Index hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            Index last_chain_pt_id = hub_chains.find(hub_id)->second.back();
            Real lastdist = my_dists[i];
            Real curdist = distance(points[i], distributed? last_chain_pt_map.find(last_chain_pt_id)->second : points[last_chain_pt_id]);

            if (curdist <= lastdist)
            {
                my_dists[i] = curdist;
                my_hub_pt_ids[i] = last_chain_pt_id;
            }
        }
    }
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::init(Real cutoff)
{
    initialize_root_hub();

    while (!hub_chains.empty())
    {
        compute_farthest_hub_pts();
        update_hub_chains();
        process_leaf_chains();
        process_split_chains();
        update_dists_and_pointers();
    }
}

template <class PointTraits, class Index>
void TreeIndex<PointTraits, Index>::radii_query(const Point& query, Real radius, IndexSet& ids) const
{
    ids.clear();
    IndexVector stack = {0};
    auto distance = PointTraits::distance();

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();

        if (distance(query, points[pt[u]]) <= radius)
            ids.insert(pt[u]);

        for (Index v : children[u])
            if (distance(query, points[pt[v]]) <= radius + max_radius * level_radius(level[v]))
                stack.push_back(v);
    }
}

template <class PointTraits, class Index>
Index TreeIndex<PointTraits, Index>::build_rgraph(Real radius, IndexSetVector& rgraph) const
{
    auto distance = PointTraits::distance();
    auto comm = points.getcomm();
    int myrank = comm.rank();
    int nprocs = comm.size();
    Index mysize = points.getmysize();

    IndexSet ids;
    Index n_edges = 0;
    rgraph.clear(); rgraph.reserve(mysize);

    assert(!comm.is_distributed()); /* TODO: temporary */

    for (Index i = 0; i < mysize; ++i)
    {
        ids.clear();
        Point p = points[i];
        radii_query(p, radius, ids);

        rgraph.push_back(ids);
        n_edges += ids.size();
    }

    return n_edges;
}
