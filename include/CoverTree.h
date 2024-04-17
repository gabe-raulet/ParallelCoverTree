#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include "Point.h"
#include <vector>

using namespace std;

class CoverTree
{
public:
    CoverTree(const vector<Point>& points, double base);

    int64_t num_points() const { return points.size(); }
    int64_t num_vertices() const { return pt.size(); }
    int64_t num_levels() const { return *max_element(level.begin(), level.end()) + 1; }

    Point get_point(int64_t point_id) const { return points[point_id]; }
    Point get_vertex_point(int64_t vertex_id) const { return points[pt[vertex_id]]; }

    vector<int64_t> radii_query(const Point& query, double radius) const;
    vector<vector<int64_t>> build_epsilon_graph(double radius) const;

    void write_gml(const char *filename) const;

private:
    double max_radius, base;
    vector<Point> points;

    vector<int64_t> pt, level;
    vector<vector<int64_t>> children;

    void build_tree();
    void set_max_radius();

    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double vertex_ball_radius(int64_t vertex_id) const;

    void update_cell(vector<double>& dists, vector<int64_t>& closest, int64_t farthest, const vector<int64_t>& descendants) const;
    vector<vector<int64_t>> compute_child_partitions(int64_t parent_id, const vector<int64_t>& descendants) const;
};

#endif
