#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include <stdint.h>
#include <math.h>
#include <vector>
#include <memory>

typedef int64_t index_t;

class CoverTree
{
public:
    CoverTree() : max_radius(0), base(2), pointmem(NULL), npoints(0), d(0) {}
    CoverTree(const float *p, index_t n, int d, double base = 2.);

    index_t num_vertices() const { return pt.size(); }
    index_t num_points() const { return npoints; }
    index_t num_levels() const { return levelset.size(); }
    int getdim() const { return d; }

    index_t get_point_id(index_t vertex_id) const { return pt[vertex_id]; }
    index_t get_vertex_level(index_t vertex_id) const { return level[vertex_id]; }
    std::vector<index_t> get_child_ids(index_t parent_id) const { return children[parent_id]; }
    std::vector<index_t> radii_query(const float *query, double radius) const;

    bool is_full() const;
    bool is_nested() const;
    bool is_covering() const;
    void print_info() const;

    void write_to_file(const char *fname) const;
    void read_from_file(const char *fname);

private:
    double max_radius, base;
    std::unique_ptr<float[]> pointmem;
    index_t npoints;
    int d;

    std::vector<index_t> pt;
    std::vector<index_t> level;
    std::vector<std::vector<index_t>> children;
    std::vector<std::vector<index_t>> levelset;

    void build_tree();
    void set_max_radius();
    index_t add_vertex(index_t point_id, index_t parent_id);

    double point_dist(index_t id1, index_t id2) const;
    double vertex_ball_radius(index_t vertex_id) const;
};

#endif
