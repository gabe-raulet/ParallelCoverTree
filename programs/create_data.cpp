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
    int nthreads = 1;
    double minval = -1.0;
    double maxval = 1.0;
    const char *output;

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -o STR    output filename [required]\n");
        fprintf(stderr, "         -d INT    point dimension [default: %d]\n", d);
        fprintf(stderr, "         -s INT    seed [default: random]\n");
        fprintf(stderr, "         -t INT    number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "         -L FLOAT  min element value [default: %.2f]\n", minval);
        fprintf(stderr, "         -U FLOAT  max element value [default: %.2f]\n", maxval);
        fprintf(stderr, "         -h        help message\n");
        return 0;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    output = read_string_arg(argc, argv, "-o", NULL);
    d = read_int_arg(argc, argv, "-d", &d);
    seed = read_int_arg(argc, argv, "-s", &seed);
    minval = read_double_arg(argc, argv, "-L", &minval);
    maxval = read_double_arg(argc, argv, "-U", &maxval);
    nthreads = read_int_arg(argc, argv, "-t", &nthreads);

    assert(minval <= maxval);

    double t;

    t = -omp_get_wtime();
    PointStore points = PointStore::pgenerate_random_points(n, d, seed, minval, maxval, nthreads);
    t += omp_get_wtime();

    fprintf(stderr, "Generated %lld random %dd points in %.4f seconds (%.2f points per millisecond)\n", n, d, t, ((n+0.0)/(t*1000)));

    t = -omp_get_wtime();
    points.write_points(output);
    t += omp_get_wtime();

    fprintf(stderr, "Wrote %lld %dd points to file '%s' in %.4f seconds (%.2f points per millisecond)\n", n, d, output, t, ((n+0.0)/(t*1000)));

    return 0;
}

