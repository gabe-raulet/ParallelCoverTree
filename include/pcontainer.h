#ifndef POINT_CONTAINER_H_
#define POINT_CONTAINER_H_

#include "mpi_env.h"

using namespace std;

template <class PointTraits, class Index>
class PointContainer
{
    public:

        using Real = typename PointTraits::Real;
        using Point = typename PointTraits::Point;
        using Distance = typename PointTraits::Distance;

        using Comm = MPIEnv::Comm;
        using PointVector = vector<Point>;
        using IndexVector = vector<Index>;

        PointContainer();
        PointContainer(const Comm& comm);
        PointContainer(const PointContainer& rhs);

        void swap(PointContainer &rhs);
        PointContainer& operator=(const PointContainer& rhs);

        void assign(const PointVector& mybuf);

        Index getmysize() const { return mysize; }
        Index gettotsize() const { return totsize; }
        Index getmyoffset() const { return myoffset; }

        bool isowner(Index globalid) const;
        Point bcast(Index globalid) const;

        int gather(PointVector& points, int root) const;
        int gather(const IndexVector& globalids, PointVector& points, int root) const;

        int allgather(PointVector& points) const;
        int allgather(const IndexVector& globalids, PointVector& points) const;
        int allgather(const IndexVector& globalids, PointVector& points, IndexVector& gatherids) const;

        void rebalance(Index chunksize);

        Point getlocal(Index localid) const { return mypoints[localid]; }
        Point operator[](Index localid) const { return mypoints[localid]; }

        Comm getcomm() const { return comm; }

    private:

        PointVector mypoints;
        Index mysize, totsize, myoffset;
        Comm comm;

        void recompute_index_info();
};

#include "pcontainer.hpp"

#endif
