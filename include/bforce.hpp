template <class PointTraits, class Real, class Index>
BForce<PointTraits, Real, Index>::BForce(MPI_Comm _comm)
    : comm(_comm),
      mysize(0),
      totsize(0),
      myoffset(0),
      distributed(false) {}

template <class PointTraits, class Real, class Index>
void BForce<PointTraits, Real, Index>::build(const PointVector& mypoints, Real cutoff)
{
    this->mypoints.assign(mypoints.begin(), mypoints.end());
    this->cutoff = cutoff;
    init();
}

template <class PointTraits, class Real, class Index>
void BForce<PointTraits, Real, Index>::init()
{
    int myrank = comm.getrank();
    int nprocs = comm.getsize();
    Distance distance = PointTraits::distance();

    distributed = !(nprocs == 1 || comm.is_self_comm());

    myoffset = 0;
    totsize = mysize = mypoints.size();
    local_cutoff_neighs.reserve(mysize);

    Index i, j;

    for (i = 0; i < mysize; ++i)
    {
        local_cutoff_neighs.emplace_back();
        IndexVector& neighs = local_cutoff_neighs.back();

        for (j = 0; j < mysize; ++j)
            if (distance(mypoints[i], mypoints[j]) <= cutoff)
                neighs.push_back(j);
    }

    if (distributed)
    {
        dist_cutoff_neighs.resize(mysize);

        comm.exscan(mysize, myoffset, MPI_SUM);
        comm.allreduce(totsize, MPI_SUM);

        PointVector allpoints;
        comm.allgatherv(mypoints, allpoints);

        i = 0;
        while (i < totsize)
        {
            if (myoffset <= i && i < myoffset + mysize)
            {
                i = myoffset + mysize;
                continue;
            }

            for (j = 0; j < mysize; ++j)
                if (distance(allpoints[i], mypoints[j]) <= cutoff)
                    dist_cutoff_neighs[j].insert({i, allpoints[i]});

            i++;
        }
    }
}

template <class PointTraits, class Real, class Index>
Index BForce<PointTraits, Real, Index>::build_rgraph(Real radius, IndexSetVector& rgraph) const
{
    int myrank = comm.getrank();
    int nprocs = comm.getsize();
    Distance distance = PointTraits::distance();

    IndexSet ids;
    Index n_edges = 0;
    rgraph.clear(); rgraph.reserve(mysize);

    for (Index i = 0; i < mysize; ++i)
    {
        ids.clear();
        Point p = mypoints[i];
        const IndexVector& local_neighs = local_cutoff_neighs[i];

        for (Index j : local_neighs)
            if (distance(p, mypoints[j]) <= radius)
                ids.insert(j + myoffset);

        if (distributed)
            for (const auto& [id, q] : dist_cutoff_neighs[i])
                if (distance(p, q) <= radius)
                    ids.insert(id);

        rgraph.push_back(ids);
        n_edges += ids.size();
    }

    return n_edges;
}

template <class PointTraits, class Real, class Index>
void BForce<PointTraits, Real, Index>::swap(BForce& other)
{
    ::swap(mypoints, other.mypoints);
    ::swap(cutoff, other.cutoff);
    ::swap(comm, other.comm);
    ::swap(dist_cutoff_neighs, other.dist_cutoff_neighs);
    ::swap(local_cutoff_neighs, other.local_cutoff_neighs);
    ::swap(mysize, other.mysize);
    ::swap(totsize, other.totsize);
    ::swap(myoffset, other.myoffset);
    ::swap(distributed, other.distributed);
}
