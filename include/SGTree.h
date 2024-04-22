#ifndef SG_TREE_H_
#define SG_TREE_H_

#include "Point.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>

using namespace std;

class SGTree
{
public:
    SGTree(const vector<Point>& points, double base);

    void build_tree();

    int64_t num_points() const;
    int64_t num_vertices() const;
    int64_t num_levels() const;

    Point get_point(int64_t point_id) const;
    Point get_vertex_point(int64_t vertex_id) const;

private:
    double max_radius, base;
    vector<Point> points;
    vector<int64_t> pt, level;
    vector<vector<int64_t>> children;

    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double vertex_ball_radius(int64_t vertex_id) const;

    /*
     * Auxiliary tree construction data structures and functions
     */
    int64_t cur_iter;
    vector<double> dists;
    vector<int64_t> hub_vtx_ids, hub_pt_ids;
    unordered_map<int64_t, vector<int64_t>> hub_chains;
    unordered_map<int64_t, int64_t> farthest_hub_pts;
    unordered_set<int64_t> split_chains, leaf_chains;

    void initialize_root_hub();
    void compute_farthest_hub_pts();
    void update_hub_chains();
    void process_leaf_chains();
    void process_split_chains();
    void update_dists_and_pointers();
};

#endif
