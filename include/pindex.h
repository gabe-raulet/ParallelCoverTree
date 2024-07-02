#ifndef POINT_INDEX_H_
#define POINT_INDEX_H_

#include "mpi_env.h"
#include "pcontainer.h"
#include <unordered_set>
#include <unordered_map>

using namespace std;

template <class PointTraits, class Index>
class PointIndex
{
    public:

        using Real = typename PointTraits::Real;
        using Point = typename PointTraits::Point;
        using Distance = typename PointTraits::Distance;

        using Comm = MPIEnv::Comm;
        using PContainer = PointContainer<PointTraits, Index>;

        using PointVector = vector<Point>;
        using IndexVector = vector<Index>;
        using IndexVectorVector = vector<IndexVector>;
        using IndexSet = unordered_set<Index>;
        using PointMap = unordered_map<Index, Point>;
        using IndexSetVector = vector<IndexSet>;
        using PointMapVector = vector<PointMap>;

        PointIndex();
        PointIndex(const Comm& comm);
        PointIndex(const PointIndex& rhs);

        void swap(PointIndex& rhs);
        PointIndex& operator=(const PointIndex& rhs);

        void build(const PointVector& mypoints, Real cutoff);

        Index getmysize() const { return points.getmysize(); }
        Index gettotsize() const { return points.gettotsize(); }
        Index getmyoffset() const { return points.getmyoffset(); }

        Index build_rgraph(Real radius, IndexSetVector& rgraph) const;

        Comm getcomm() const { return points.getcomm(); }

    private:

        PContainer points;
        IndexVectorVector local_neighs;
        PointMapVector dist_neighs;

        void init(Real cutoff);
};

#include "pindex.hpp"

#endif
