#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include "Point.h"
#include <vector>
#include <cstdint>
#include <cstddef>

class CoverTree
{
public:
    CoverTree(const std::vector<Point>& points) { build_tree(); }

    size_t num_vertices() const { return pt.size(); }
    size_t num_points() const { return points.size(); }
    size_t num_levels() const { return levelset.size(); }

    size_t radii_query(const Point& query, double radius, std::vector<int64_t>& ids) const; /* TODO */

    int64_t get_vertex_level(int64_t vertex_id) const { return level[vertex_id]; }
    int64_t get_point_id(int64_t vertex_id) const { return pt[vertex_id]; }
    int64_t get_parent_id(int64_t vertex_id) const { return parent[vertex_id]; }

    size_t get_child_ids(int64_t vertex_id, std::vector<int64_t>& child_ids) const;
    size_t get_level_ids(int64_t vertex_level, std::vector<int64_t>& level_ids) const;

private:
    std::vector<Point> points;
    double max_radius;

    std::vector<int64_t> pt;
    std::vector<int64_t> parent;
    std::vector<int64_t> level;

    std::vector<std::vector<int64_t>> children;
    std::vector<std::vector<int64_t>> levelset;

    void build_tree(); /* TODO */
    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double point_dist(int64_t point_id1, int64_t point_id2) const;
    double vertex_ball_radius(int64_t vertex_id) const;
};

#endif
