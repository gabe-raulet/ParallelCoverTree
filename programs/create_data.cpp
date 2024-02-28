#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <cassert>
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
    int64_t n;
    int d = 2;
    int seed = -1;
    const char *output;

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -o STR    output filename [required]\n");
        fprintf(stderr, "         -d INT    point dimension [default: %d]\n", d);
        fprintf(stderr, "         -s INT    seed [default: random]\n");
        fprintf(stderr, "         -h        help message\n");
        return 0;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    output = read_string_arg(argc, argv, "-o", NULL);
    d = read_int_arg(argc, argv, "-d", &d);
    seed = read_int_arg(argc, argv, "-s", &seed);

    double t;

    t = -omp_get_wtime();
    PointStore points = PointStore::generate_random_points(n, d, seed, -1.0, 1.0);
    //std::vector<Point> points = Point::random_points(n, d, seed);
    t += omp_get_wtime();

    fprintf(stderr, "Generated %lld random %dd points in %.4f seconds (%.2f points per millisecond)\n", n, d, t, ((n+0.0)/(t*1000)));

    t = -omp_get_wtime();
    points.write_points(output);
    //Point::write_points(points, output);
    t += omp_get_wtime();

    fprintf(stderr, "Wrote %lld %dd points to file '%s' in %.4f seconds (%.2f points per millisecond)\n", n, d, output, t, ((n+0.0)/(t*1000)));

    return 0;
}

