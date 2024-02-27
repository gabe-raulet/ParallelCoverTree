#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <omp.h>

typedef std::vector<std::vector<int64_t>> nng_t;

bool nngs_are_equal(const nng_t& a, const nng_t& b);

int main(int argc, char *argv[])
{
    int max_nthrds = omp_get_max_threads();

    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <points.bin> <radius>\n", argv[0]);
        return 0;
    }

    const char *points_fname = argv[1];
    double radius = atof(argv[2]);

    double telapsed, t_serial_build, t_serial_nng;

    std::vector<Point> points;
    nng_t truth, nng;

    Point::read_points(points, points_fname);

    omp_set_num_threads(1);

    t_serial_build = -omp_get_wtime();
    CoverTree ct = CoverTree::build_recursive(points);
    t_serial_build += omp_get_wtime();

    printf("Successfully built cover tree serially in %.4f seconds\n", t_serial_build);

    t_serial_nng = ct.get_neighborhood_graph(radius, truth);

    printf("Built neighborhood graph serially in %.4f seconds\n\n", t_serial_nng);

    for (int nthrds = 2; nthrds <= max_nthrds; nthrds *= 2)
    {
        omp_set_num_threads(nthrds);

        telapsed = -omp_get_wtime();
        CoverTree tree = CoverTree::build_recursive(points);
        telapsed += omp_get_wtime();

        printf("Successfully built cover tree with %d threads in %.4f seconds (%.2f speedup)\n", nthrds, telapsed, t_serial_build/telapsed);

        telapsed = ct.get_neighborhood_graph(radius, nng);
        assert(nngs_are_equal(truth, nng));

        printf("Successfully built neighborhood graph with %d threads in %.4f seconds (%.2f speedup)\n\n", nthrds, telapsed, t_serial_nng/telapsed);
    }

    return 0;
}

bool nngs_are_equal(const nng_t& a, const nng_t& b)
{
    uint64_t n = a.size();

    if (n != b.size())
        return false;

    auto cmp1 = a.cbegin();
    auto cmp2 = b.cbegin();

    for (uint64_t i = 0; i < n; ++i)
        if (*cmp1++ != *cmp2++)
            return false;

    return true;
}
