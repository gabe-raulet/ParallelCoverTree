#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include "Point.h"
#include <vector>

using namespace std;

class CoverTree
{
public:
    CoverTree(const vector<Point>& points, double base);

    void build_tree_hub_loop();
    void build_tree_point_loop();

    int64_t num_points() const;
    int64_t num_vertices() const;
    int64_t num_levels() const;

    Point get_point(int64_t point_id) const;
    Point get_vertex_point(int64_t vertex_id) const;

    vector<int64_t> radii_query(const Point& query, double radius) const;
    vector<vector<int64_t>> build_epsilon_graph(double radius) const;

    void write_gml(const char *filename) const;

private:
    double max_radius, base;
    vector<Point> points;

    vector<int64_t> pt, level;
    vector<vector<int64_t>> children;

    void set_max_radius();

    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double vertex_ball_radius(int64_t vertex_id) const;

    void update_cell(vector<double>& dists, vector<int64_t>& closest, int64_t farthest, const vector<int64_t>& descendants) const;
    vector<vector<int64_t>> compute_child_partitions(int64_t parent_id, const vector<int64_t>& descendants) const;
};

#endif
