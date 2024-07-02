template <class PointTraits, class Index>
PointIndex<PointTraits, Index>::PointIndex()
    : points(Comm::comm_self()) {}

template <class PointTraits, class Index>
PointIndex<PointTraits, Index>::PointIndex(const Comm& comm)
    : points(comm) {}

template <class PointTraits, class Index>
PointIndex<PointTraits, Index>::PointIndex(const PointIndex& rhs)
    : points(rhs.points) {}

template <class PointTraits, class Index>
void PointIndex<PointTraits, Index>::swap(PointIndex& rhs)
{
    points.swap(rhs.points);
}

template <class PointTraits, class Index>
PointIndex<PointTraits, Index>& PointIndex<PointTraits, Index>::operator=(const PointIndex& rhs)
{
    PointIndex tmp(rhs);
    tmp.swap(*this);
    return *this;
}

template <class PointTraits, class Index>
void PointIndex<PointTraits, Index>::build(const PointVector& mypoints, Real cutoff)
{
    points.assign(mypoints);
    init(cutoff);
}

template <class PointTraits, class Index>
void PointIndex<PointTraits, Index>::init(Real cutoff)
{
    auto distance = PointTraits::distance();
    auto comm = points.getcomm();
    Index mysize = getmysize();
    Index myoffset = getmyoffset();
    Index totsize = gettotsize();

    local_neighs.clear();
    local_neighs.reserve(mysize);

    Index i, j;

    for (i = 0; i < mysize; ++i)
    {
        local_neighs.emplace_back();
        auto& neighs = local_neighs.back();

        for (j = 0; j < mysize; ++j)
            if (distance(points[i], points[j]) <= cutoff)
                neighs.push_back(j);
    }

    if (comm.is_distributed())
    {
        dist_neighs.resize(mysize);
        comm.exscan(mysize, myoffset, MPI_SUM, static_cast<Index>(0));

        PointVector allpoints;
        points.allgather(allpoints);

        i = 0;
        while (i < totsize)
        {
            if (myoffset <= i && i < myoffset + mysize)
            {
                i = myoffset + mysize;
                continue;
            }

            for (j = 0; j < mysize; ++j)
                if (distance(allpoints[i], points[j]) <= cutoff)
                    dist_neighs[j].insert({i, allpoints[i]});

            i++;
        }
    }
}

template <class PointTraits, class Index>
Index PointIndex<PointTraits, Index>::build_rgraph(Real radius, IndexSetVector& rgraph) const
{
    auto comm = points.getcomm();
    int myrank = comm.rank();
    int nprocs = comm.size();
    auto distance = PointTraits::distance();

    IndexSet ids;
    Index n_edges = 0;
    Index mysize = getmysize();
    Index myoffset = getmyoffset();
    rgraph.clear(); rgraph.reserve(mysize);

    for (Index i = 0; i < mysize; ++i)
    {
        ids.clear();
        auto p = points[i];
        const auto& neighs = local_neighs[i];

        for (Index j : neighs)
            if (distance(p, points[j]) <= radius)
                ids.insert(myoffset + j);

        if (comm.is_distributed())
            for (const auto& [id, q] : dist_neighs[i])
                if (distance(p, q) <= radius)
                    ids.insert(id);

        rgraph.push_back(ids);
        n_edges += ids.size();
    }

    return n_edges;
}
