#ifndef DIST_COVER_TREE_H_
#define DIST_COVER_TREE_H_

#include "Point.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <mpi.h>

using namespace std;

class DistCoverTree
{
public:
    DistCoverTree(const vector<Point>& mypoints, double base, MPI_Comm comm);

    void build_tree();

    int64_t my_num_points() const { return mysize; }
    int64_t tot_num_points() const { return totsize; }
    int64_t my_points_offset() const { return myoffset; }

    Point get_my_point(int64_t point_id) const;
    Point get_my_vertex_point(int64_t vertex_id) const;

    int64_t num_vertices() const { return 0; }
    int64_t num_levels() const {return 0; }

private:
    double max_radius, base;
    vector<int64_t> pt, level;
    vector<vector<int64_t>> children;

    vector<Point> mypoints;
    int64_t mysize, totsize, myoffset;
    MPI_Comm comm;

    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double vertex_ball_radius(int64_t vertex_id) const;

    vector<double> my_dists;
    vector<int64_t> my_hub_vtx_ids, my_hub_pt_ids;
    unordered_map<int64_t, vector<int64_t>> hub_chains;
    unordered_map<int64_t, pair<int64_t, double>> farthest_hub_pts;
    unordered_set<int64_t> leaf_chains;
    vector<int64_t> split_chains;

    void initialize_root_hub();
    void compute_farthest_hub_pts();
    void update_hub_chains();
    void process_leaf_chains();
    void process_split_chains();
    void update_dists_and_pointers();
};

#endif
