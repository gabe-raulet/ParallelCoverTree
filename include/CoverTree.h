#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include "Point.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;

class CoverTree
{
public:
    CoverTree(const vector<Point>& points, double base);

    void build_tree(bool verbose = false);

    int64_t num_points() const;
    int64_t num_vertices() const;
    int64_t num_levels() const;

    Point get_point(int64_t point_id) const;
    Point get_vertex_point(int64_t vertex_id) const;

    vector<int64_t> radii_query(const Point& query, double radius) const;
    vector<vector<int64_t>> build_epsilon_graph(double radius) const;

    void print_timing_results() const;

private:
    double max_radius, base;
    vector<Point> points;
    vector<int64_t> pt, level;
    vector<vector<int64_t>> children;
    int64_t nlevels, niters, num_active_pts;

    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double vertex_ball_radius(int64_t vertex_id) const;

    /*
     * Auxiliary tree construction data structures and functions
     */
    int64_t cur_iter;
    vector<double> dists;
    vector<int64_t> hub_vtx_ids, hub_pt_ids;
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
};

#endif
