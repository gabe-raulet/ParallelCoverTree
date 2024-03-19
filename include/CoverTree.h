#ifndef COVER_TREE_H_
#define COVER_TREE_H_

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <memory>

typedef int64_t index_t;

class CoverTree
{
public:
    CoverTree(const float *p, index_t n, int d, double base = 2.);

    index_t num_points() const { return npoints; }
    index_t num_vertices() const { return pt.size(); }
    int getdim() const { return d; }
    const float* getdata() const { return pointmem.get(); }

    std::vector<index_t> radii_query(const float *query, double radius) const;

    bool is_full() const;
    bool is_nested() const;
    bool is_covering() const;
    void print_info(FILE *f) const;

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

    void build_tree();
    void set_max_radius();
    index_t add_vertex(index_t point_id, index_t parent_id);

    double point_dist(index_t id1, index_t id2) const;
    double vertex_ball_radius(index_t vertex_id) const;

    std::vector<std::vector<index_t>> get_level_set() const;

    std::vector<std::tuple<index_t, std::vector<index_t>>>
    compute_child_points(index_t parent_id, const std::vector<index_t>& descendants) const;
};

#endif
