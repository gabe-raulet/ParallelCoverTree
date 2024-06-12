template <class PointTraits, class Real, class Index>
CoverTree<PointTraits, Real, Index>::CoverTree(MPI_Comm _comm)
    : comm(_comm),
      mysize(0),
      totsize(0),
      myoffset(0),
      distributed(false) {}

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::build(const PointVector& mypoints, Real cutoff, Real base)
{
    this->mypoints.assign(mypoints.begin(), mypoints.end());
    this->cutoff = cutoff;
    this->base = base;

    int myrank = comm.getrank();
    int nprocs = comm.getsize();

    distributed = !(nprocs == 1 || comm.is_self_comm());

    myoffset = 0;
    totsize = mysize = mypoints.size();

    if (distributed)
    {
        comm.exscan(mysize, myoffset, MPI_SUM);
        comm.allreduce(totsize, MPI_SUM);
    }

    init();
}

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::init()
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

template <class PointTraits, class Real, class Index>
Index CoverTree<PointTraits, Real, Index>::add_vertex(Index point, Index parent)
{
    Index vertex_level;
    Index vertex = pt.size();

    pt.push_back(point);
    children.emplace_back();

    if (parent >= 0)
    {
        vertex_level = level[parent] + 1;
        children[parent].push_back(vertex);
    }
    else vertex_level = 0;

    nlevels = max(vertex_level+1, nlevels);
    level.push_back(vertex_level);

    return vertex;
}

template <class PointTraits, class Real, class Index>
Real CoverTree<PointTraits, Real, Index>::level_radius(Index level) const
{
    return pow(base, -1. * level);
}

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::initialize_root_hub()
{
    my_dists.resize(mysize);
    my_hub_vtx_ids.resize(mysize);
    my_hub_pt_ids.resize(mysize);

    Index root = add_vertex(0, -1);
    hub_chains.insert({root, {pt[root]}});

    Distance distance = PointTraits::distance();
    Point root_pt = mypoints.front();

    if (distributed) comm.bcast(root_pt, 0);

    for (Index i = 0; i < mysize; ++i)
    {
        my_dists[i] = distance(root_pt, mypoints[i]);
        my_hub_vtx_ids[i] = root;
        my_hub_pt_ids[i] = pt[root];
        max_radius = max(max_radius, my_dists[i]);
    }

    if (distributed) comm.allreduce(max_radius, MPI_MAX);
}

#include "mpi_argmax.h"

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::compute_farthest_hub_pts()
{
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

    if (!distributed)
    {
        farthest_hub_pts = ::move(my_argmaxes);
        return;
    }

    using ArgmaxPair = ArgmaxPair<Real, Index>;
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

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::update_hub_chains()
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

template <class PointTraits, class Real, class Index>
Index CoverTree<PointTraits, Real, Index>::batch_new_vertex(Index point, Index parent)
{
    my_new_vertex_pt_ids.push_back(point);
    my_new_vertex_hub_ids.push_back(parent);
    return pt.size() + my_new_vertex_pt_ids.size();
}

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::add_batched_vertices()
{
    assert(my_new_vertex_pt_ids.size() == my_new_vertex_hub_ids.size());

    if (!distributed)
        for (Index i = 0; i < my_new_vertex_pt_ids.size(); ++i)
            add_vertex(my_new_vertex_pt_ids[i], my_new_vertex_hub_ids[i]);
    else
    {
        IndexVector new_vertex_pt_ids, new_vertex_hub_ids;

        comm.allgatherv(my_new_vertex_pt_ids, new_vertex_pt_ids);
        comm.allgatherv(my_new_vertex_hub_ids, new_vertex_hub_ids);

        for (Index i = 0; i < new_vertex_pt_ids.size(); ++i)
            add_vertex(new_vertex_pt_ids[i], new_vertex_hub_ids[i]);

        my_new_vertex_pt_ids.clear();
        my_new_vertex_hub_ids.clear();
    }
}

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::process_leaf_chains()
{
    if (leaf_chains.empty())
        return;

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

    add_batched_vertices();
}

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::process_split_chains()
{
    using IndexMap = unordered_map<Index, Index>;

    if (split_chains.empty())
        return;

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

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::update_dists_and_pointers()
{
    Distance distance = PointTraits::distance();
    PointMap last_chain_pt_map;

    if (distributed)
    {
        PointVector my_last_chain_pts;
        IndexVector my_last_chain_pt_ids;

        for (auto it = hub_chains.begin(); it != hub_chains.end(); ++it)
        {
            Index last_chain_pt_id = it->second.back();

            if (myoffset <= last_chain_pt_id && last_chain_pt_id < myoffset + mysize)
            {
                my_last_chain_pt_ids.push_back(last_chain_pt_id);
                my_last_chain_pts.push_back(mypoints[last_chain_pt_id - myoffset]);
            }
        }

        PointVector last_chain_pts;
        IndexVector last_chain_pt_ids;

        comm.allgatherv(my_last_chain_pts, last_chain_pts);
        comm.allgatherv(my_last_chain_pt_ids, last_chain_pt_ids);

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
            Real curdist = distance(mypoints[i], distributed? last_chain_pt_map.find(last_chain_pt_id)->second : mypoints[last_chain_pt_id]);

            if (curdist <= lastdist)
            {
                my_dists[i] = curdist;
                my_hub_pt_ids[i] = last_chain_pt_id;
            }
        }
    }
}

template <class PointTraits, class Real, class Index>
void CoverTree<PointTraits, Real, Index>::radii_query(const Point& query, Real radius, IndexSet& ids) const
{
    ids.clear();
    IndexVector stack = {0};
    Distance distance = PointTraits::distance();

    while (!stack.empty())
    {
        Index u = stack.back(); stack.pop_back();

        if (distance(query, mypoints[pt[u]]) <= radius)
            ids.insert(pt[u]);

        for (Index v : children[u])
            if (distance(query, mypoints[pt[v]]) <= radius + max_radius * level_radius(level[v]))
                stack.push_back(v);
    }
}

template <class PointTraits, class Real, class Index>
Index CoverTree<PointTraits, Real, Index>::build_rgraph(Real radius, IndexSetVector& rgraph) const
{
    int myrank = comm.getrank();
    int nprocs = comm.getsize();
    Distance distance = PointTraits::distance();

    IndexSet ids;
    Index n_edges = 0;
    rgraph.clear(); rgraph.reserve(mysize);

    assert(!distributed); /* TODO: temporary */

    for (Index i = 0; i < mysize; ++i)
    {
        ids.clear();
        Point p = mypoints[i];
        radii_query(p, radius, ids);

        rgraph.push_back(ids);
        n_edges += ids.size();
    }

    return n_edges;
}
