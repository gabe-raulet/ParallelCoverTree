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
    int max_nthrds = omp_get_max_threads();

    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <points.bin>\n", argv[0]);
        return 0;
    }

    const char *points_fname = argv[1];

    double telapsed, t_serial_build;

    //std::vector<Point> points;
    PointStore points = PointStore::read_points(points_fname);

    omp_set_num_threads(1);

    t_serial_build = -omp_get_wtime();
    CoverTree ct = CoverTree::build(points);
    t_serial_build += omp_get_wtime();

    printf("Successfully built cover tree serially in %.4f seconds\n\n", t_serial_build);
    ct.print_info();

    return 0;
}
