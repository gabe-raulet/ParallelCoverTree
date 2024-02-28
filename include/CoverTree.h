#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include "Point.h"
#include <vector>
#include <cstdint>
#include <cstddef>

class CoverTree
{
public:
    CoverTree(const std::vector<Point>& points) : points(points) {}

    static CoverTree build(const std::vector<Point>& points)
    {
        CoverTree ct(points);
        ct.build_tree();
        return ct;
    }

    static CoverTree build_recursive(const std::vector<Point>& points)
    {
        CoverTree ct(points);
        ct.build_tree_recursive();
        return ct;
    }

    size_t num_vertices() const { return pt.size(); }
    size_t num_points() const { return points.size(); }
    size_t num_levels() const { return levelset.size(); }

    size_t radii_query(const Point& query, double radius, std::vector<int64_t>& ids) const;
    double get_neighborhood_graph(double radius, std::vector<std::vector<int64_t>>& nng, bool sort = true) const;

    int64_t get_vertex_level(int64_t vertex_id) const { return level[vertex_id]; }
    int64_t get_point_id(int64_t vertex_id) const { return pt[vertex_id]; }
    int64_t get_parent_id(int64_t vertex_id) const { return parent[vertex_id]; }

    size_t get_child_ids(int64_t vertex_id, std::vector<int64_t>& child_ids) const;
    size_t get_level_ids(int64_t vertex_level, std::vector<int64_t>& level_ids) const;

    const std::vector<Point>& get_points() const { return points; }

    double average_vertex_degree(int64_t level) const;
    void print_info() const;

private:
    std::vector<Point> points;
    double max_radius;

    std::vector<int64_t> pt;
    std::vector<int64_t> parent;
    std::vector<int64_t> level;

    std::vector<std::vector<int64_t>> children;
    std::vector<std::vector<int64_t>> levelset;

    void build_tree();
    void build_tree_recursive();
    void build_tree_recursive_f(int64_t parent_id, const std::vector<int64_t>& descendants);
    int64_t add_vertex(int64_t point_id, int64_t parent_id);
    double point_dist(int64_t point_id1, int64_t point_id2) const;
    double vertex_ball_radius(int64_t vertex_id) const;

    std::vector<std::tuple<int64_t, std::vector<int64_t>>>
    get_next_parents(int64_t parent_id, const std::vector<int64_t>& descendants);

    bool is_full() const;
    bool is_nested() const;
    bool is_covering() const;

    void set_max_radius();
    std::vector<std::tuple<int64_t, std::vector<int64_t>>> init_build_stack();
};

#endif
