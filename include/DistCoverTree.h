#ifndef DIST_COVER_TREE_H_
#define DIST_COVER_TREE_H_

#include "Point.h"
#include "CoverTree.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <mpi.h>

using namespace std;

class DistCoverTree
{
public:
    DistCoverTree(const vector<Point>& mypoints, double base, MPI_Comm comm);

    void build_tree(bool verbose = false);

    int64_t num_vertices() const;
    int64_t num_levels() const;

    void print_timing_results() const;

    vector<vector<int64_t>> build_epsilon_graph(double radius) const;

private:
    double max_radius, base;
    vector<int64_t> pt, level;
    unordered_map<int64_t, Point> repoints;
    vector<vector<int64_t>> children;
    unordered_map<int64_t, double> cover_map;
    int64_t niters, nlevels;

    vector<Point> mypoints;
    int64_t mysize, totsize, myoffset;
    MPI_Comm comm;

    unordered_map<int64_t, vector<int64_t>> local_ptid_map;
    unordered_map<int64_t, CoverTree> local_trees;

    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double vertex_ball_radius(int64_t vertex_id) const;

    vector<double> my_dists;
    vector<int64_t> my_hub_vtx_ids, my_hub_pt_ids;
    unordered_map<int64_t, vector<int64_t>> hub_chains;
    unordered_map<int64_t, pair<int64_t, double>> farthest_hub_pts;
    unordered_set<int64_t> leaf_chains;
    vector<int64_t> split_chains;

    void initialize_root_hub(bool verbose = false);
    void compute_farthest_hub_pts(bool verbose = false);
    void update_hub_chains(bool verbose = false);
    void process_leaf_chains(bool verbose = false);
    void process_split_chains(bool verbose = false);
    void update_dists_and_pointers(bool verbose = false);

    double overall_time;
    double initialize_root_hub_time;
    double compute_farthest_hub_pts_time;
    double update_hub_chains_time;
    double process_leaf_chains_time;
    double process_split_chains_time;
    double update_dists_and_pointers_time;

    void set_times_to_zero();

    vector<int64_t> my_new_vertex_pt_ids, my_new_vertex_hub_ids;
    int64_t batch_new_vertex(int64_t point_id, int64_t parent_id);
    void add_batched_vertices();

    unordered_map<int64_t, int64_t> get_hub_counts() const;
    pair<double, unordered_map<int64_t, int>> compute_hub_assignments(bool verbose) const;
    void collect_replicate_points(bool verbose = false);
    void build_local_trees(const unordered_map<int64_t, int>& hub_assignments, bool verbose = false);

    void dump_info() const;

    vector<int64_t> my_radii_query(const Point& query, double radius) const;
    unordered_set<int64_t> local_radii_query(const Point& query, double radius, vector<int64_t>& local_tree_ids) const;
};

#endif
