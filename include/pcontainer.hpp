template <class PointTraits, class Index>
PointContainer<PointTraits, Index>::PointContainer()
    : comm(Comm::comm_self()) {}

template <class PointTraits, class Index>
PointContainer<PointTraits, Index>::PointContainer(const Comm& comm)
    : comm(comm) {}

template <class PointTraits, class Index>
PointContainer<PointTraits, Index>::PointContainer(const PointContainer& rhs)
    : mypoints(rhs.mypoints),
      comm(rhs.comm) { recompute_index_info(); }

template <class PointTraits, class Index>
void PointContainer<PointTraits, Index>::swap(PointContainer &rhs)
{
    ::swap(mypoints, rhs.mypoints);
    ::swap(mysize, rhs.mysize);
    ::swap(totsize, rhs.totsize);
    ::swap(myoffset, rhs.myoffset);
    comm.swap(rhs.comm);
}

template <class PointTraits, class Index>
PointContainer<PointTraits, Index>&
PointContainer<PointTraits, Index>::operator=(const PointContainer& rhs)
{
    PointContainer tmp(rhs);
    tmp.swap(*this);
    return *this;
}

template <class PointTraits, class Index>
void PointContainer<PointTraits, Index>::assign(const PointVector& mybuf)
{
    mypoints.assign(mybuf.begin(), mybuf.end());
    recompute_index_info();
}

template <class PointTraits, class Index>
bool PointContainer<PointTraits, Index>::isowner(Index globalid) const
{
    return (myoffset <= globalid && globalid < myoffset + mysize);
}

template <class PointTraits, class Index>
typename PointTraits::Point
PointContainer<PointTraits, Index>::bcast(Index globalid) const
{
    /*
     * TODO: this really ought to be implemented using send/recv...
     */

    Point p;
    int root = -1;

    if (isowner(globalid))
    {
        root = comm.rank();
        p = mypoints[globalid-myoffset];
    }

    comm.allreduce(root, MPI_MAX);
    comm.bcast(p, root);

    return p;
}

template <class PointTraits, class Index>
int PointContainer<PointTraits, Index>::gather(PointVector& points, int root) const
{
    return comm.gatherv(mypoints, points, root);
}

template <class PointTraits, class Index>
int PointContainer<PointTraits, Index>::gather(const IndexVector& globalids, PointVector& points, int root) const
{
    PointVector mybuf;

    for (Index globalid : globalids)
        if (isowner(globalid))
            mybuf.push_back(mypoints[globalid-myoffset]);

    return comm.gatherv(mybuf, points, root);
}

template <class PointTraits, class Index>
int PointContainer<PointTraits, Index>::allgather(PointVector& points) const
{
    return comm.allgatherv(mypoints, points);
}

template <class PointTraits, class Index>
int PointContainer<PointTraits, Index>::allgather(const IndexVector& globalids, PointVector& points) const
{
    PointVector mybuf;

    for (Index globalid : globalids)
        if (isowner(globalid))
            mybuf.push_back(mypoints[globalid-myoffset]);

    return comm.allgatherv(mybuf, points);
}

template <class PointTraits, class Index>
int PointContainer<PointTraits, Index>::allgather(const IndexVector& globalids, PointVector& points, IndexVector& gatherids) const
{
    PointVector mybuf;
    IndexVector mygatherids;

    for (Index globalid : globalids)
        if (isowner(globalid))
        {
            mybuf.push_back(mypoints[globalid-myoffset]);
            mygatherids.push_back(globalid);
        }

    int err = comm.allgatherv(mybuf, points);
    if (err != MPI_SUCCESS) return err;

    return comm.allgatherv(mygatherids, gatherids);
}

template <class PointTraits, class Index>
void PointContainer<PointTraits, Index>::rebalance(Index chunksize)
{
    /* TODO: */
}

template <class PointTraits, class Index>
void PointContainer<PointTraits, Index>::recompute_index_info()
{
    mysize = totsize = mypoints.size();
    comm.allreduce(totsize, MPI_SUM);
    comm.exscan(mysize, myoffset, MPI_SUM, static_cast<Index>(0));
}
