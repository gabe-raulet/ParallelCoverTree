#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "mpi_comm.h"

using namespace std;

template <class PointTraits, class Real = double, class Index = int64_t>
class CoverTree
{
    public:

        using Point = typename PointTraits::Point;
        using Distance = typename PointTraits::Distance;

        using PointVector = vector<Point>;
        using IndexVector = vector<Index>;
        using RealVector = vector<Real>;
        using IndexVectorVector = vector<IndexVector>;

        using DistPair = pair<Index, Real>;

        using IndexVectorMap = unordered_map<Index, IndexVector>;
        using PointMap = unordered_map<Index, Point>;
        using IndexSet = unordered_set<Index>;
        using DistPairMap = unordered_map<Index, DistPair>;
        using IndexSetVector = vector<IndexSet>;

        CoverTree(MPI_Comm comm);
        void build(const PointVector& mypoints, Real cutoff, Real base);

        MPI_Comm getcomm() const { return comm.getcomm(); }

        Index getmysize() const { return mysize; }
        Index gettotsize() const { return totsize; }
        Index getmyoffset() const { return myoffset; }

        void radii_query(const Point& query, Real radius, IndexSet& ids) const;
        Index build_rgraph(Real radius, IndexSetVector& rgraph) const;

    private:

        Real max_radius, base;
        IndexVector pt, level;
        IndexVectorVector children;
        Index nlevels;

        PointVector mypoints;
        Index mysize, myoffset, totsize;
        Real cutoff;
        MPICommunicator comm;
        bool distributed;

        Index add_vertex(Index point, Index parent);
        Real level_radius(Index level) const;

        RealVector my_dists;
        IndexVector my_hub_vtx_ids, my_hub_pt_ids;
        IndexVectorMap hub_chains;
        DistPairMap farthest_hub_pts;
        IndexSet leaf_chains;
        IndexVector split_chains;

        void initialize_root_hub();
        void compute_farthest_hub_pts();
        void update_hub_chains();
        void process_leaf_chains();
        void process_split_chains();
        void update_dists_and_pointers();

        IndexVector my_new_vertex_pt_ids, my_new_vertex_hub_ids;
        Index batch_new_vertex(Index point, Index parent);
        void add_batched_vertices();

        void init();
};

#include "cover_tree.hpp"

#endif
