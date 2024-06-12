#ifndef BRUTE_FORCE_H_
#define BRUTE_FORCE_H_

#include <vector>
#include <unordered_set>
#include <type_traits>
#include "mpi_comm.h"

using namespace std;

template <class PointTraits, class Real = double, class Index = int64_t>
class BForce
{
    public:

        static_assert(is_same_v<Real, double> && is_same_v<Index, int64_t>); // temporary

        using Point    = typename PointTraits::Point;
        using Distance = typename PointTraits::Distance;

        using IndexSet = unordered_set<Index>;
        using PointMap = unordered_map<Index, Point>;
        using IndexSetVector = vector<IndexSet>;

        using PointVector       = vector<Point>;
        using IndexVector       = vector<Index>;
        using IndexVectorVector = vector<IndexVector>;
        using PointMapVector    = vector<PointMap>;

        BForce(MPI_Comm comm);
        void build(const PointVector& mypoints, Real cutoff);

        Index getmysize() const { return mysize; }
        Index gettotsize() const { return totsize; }
        Index getmyoffset() const { return myoffset; }

        Index build_rgraph(Real radius, IndexSetVector& rgraph) const;

        void swap(BForce& other);

        MPI_Comm getcomm() const { return comm.getcomm(); }

    private:

        PointVector mypoints;
        Real cutoff;
        MPICommunicator comm;
        PointMapVector dist_cutoff_neighs;
        IndexVectorVector local_cutoff_neighs;
        Index mysize, totsize, myoffset;
        bool distributed;

        void init();
};

#include "bforce.hpp"

#endif
