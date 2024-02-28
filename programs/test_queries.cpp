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

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <points.bin> <radius>\n", argv[0]);
        return 0;
    }

    const char *points_fname = argv[1];
    double radius = atof(argv[2]);

    double t, t2;
    //std::vector<Point> points;
    std::vector<std::vector<int64_t>> truth, nng;

    //Point::read_points(points, points_fname);
    PointStore points = PointStore::read_points(points_fname);
    CoverTree ct = CoverTree::build(points);

    omp_set_num_threads(1);
    t = ct.get_neighborhood_graph(radius, truth);

    for (int nthrds = 2; nthrds <= 64; nthrds *= 2)
    {
        omp_set_num_threads(nthrds);
        t2 = ct.get_neighborhood_graph(radius, nng);

        auto cmp1 = truth.begin();
        auto cmp2 = nng.begin();

        for (uint64_t i = 0; i < ct.num_points(); ++i)
            assert(*cmp1++ == *cmp2++);

        fprintf(stderr, "Successfully built neighborhood graph in %.5f seconds using %d threads (%.5f speedup)\n", t2, nthrds, t/t2);
    }

    return 0;
}

