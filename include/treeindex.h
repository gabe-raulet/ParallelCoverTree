#ifndef COVER_TREE_INDEX_H_
#define COVER_TREE_INDEX_H_

#include "mpi_env.h"
#include "pcontainer.h"
#include <unordered_set>
#include <unordered_map>

using namespace std;

template <class PointTraits, class Index>
class TreeIndex
{
    public:

        using Real = typename PointTraits::Real;
        using Point = typename PointTraits::Point;
        using Distance = typename PointTraits::Distance;

        using Comm = MPIEnv::Comm;
        using PContainer = PointContainer<PointTraits, Index>;

        using DistPair = pair<Index, Real>;
        using PointVector = vector<Point>;
        using IndexVector = vector<Index>;
        using RealVector = vector<Real>;
        using IndexVectorVector = vector<IndexVector>;
        using IndexSet = unordered_set<Index>;
        using PointMap = unordered_map<Index, Point>;
        using IndexVectorMap = unordered_map<Index, IndexVector>;
        using DistPairMap = unordered_map<Index, DistPair>;
        using IndexSetVector = vector<IndexSet>;
        using PointMapVector = vector<PointMap>;

        TreeIndex();
        TreeIndex(const Comm& comm);
        TreeIndex(const TreeIndex& rhs);

        void swap(TreeIndex& rhs);
        TreeIndex& operator=(const TreeIndex& rhs);

        void build(const PointVector& mypoints, Real cutoff, Real base = 2.0);

        Index build_rgraph(Real radius, IndexSetVector& rgraph) const;
        void radii_query(const Point& query, Real radius, IndexSet& ids) const;

        Comm getcomm() const { return points.getcomm(); }

    private:

        Real base; /* cover tree ball radii decay factor */
        Real max_radius; /* distance between point 0 and farthest point */

        IndexVector pt; /* maps vertex ids to point ids */
        IndexVector level; /* maps vertex ids to vertex level */
        IndexVectorVector children; /* maps vertex ids to chilren ids */
        Index nlevels; /* current number levels in tree */

        RealVector my_dists;
        IndexVector my_hub_vtx_ids;
        IndexVector my_hub_pt_ids;
        IndexVectorMap hub_chains;
        DistPairMap farthest_hub_pts;
        IndexSet leaf_chains;
        IndexVector split_chains;

        PContainer points; /* distributed point container */

        IndexVector my_new_vtx_pt_ids;
        IndexVector my_new_vtx_hub_ids;

        Index add_vertex(Index point, Index parent);
        Index batch_new_vertex(Index point, Index parent);
        void add_batched_vertices();

        Real level_radius(Index level) const { return ::pow(base, -1. * level); }

        void init(Real cutoff);
        void initialize_root_hub();
        void compute_farthest_hub_pts();
        void update_hub_chains();
        void process_leaf_chains();
        void process_split_chains();
        void update_dists_and_pointers();
};

#include "treeindex.hpp"

#endif
