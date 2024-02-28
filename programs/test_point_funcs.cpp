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

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -d INT    point dimension [default: %d]\n", d);
        fprintf(stderr, "         -s INT    seed [default: random]\n");
        fprintf(stderr, "         -t INT    number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "         -h        help message\n");
        return 0;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    d = read_int_arg(argc, argv, "-d", &d);
    seed = read_int_arg(argc, argv, "-s", &seed);
    nthreads = read_int_arg(argc, argv, "-t", &nthreads);

    fprintf(stderr, "** Using %d threads **\n", std::min(nthreads, omp_get_max_threads()));

    double t;

    t = -omp_get_wtime();
    PointStore A = PointStore::pgenerate_random_points(n, d, seed, -1.0, 1.0, nthreads);
    t += omp_get_wtime();

    fprintf(stderr, "Generated %lld random %dd points in %.4f seconds (%.2f points per millisecond)\n", n, d, t, ((n+0.0)/(t*1000)));

    t = -omp_get_wtime();
    A.write_points("points.bin");
    t += omp_get_wtime();

    fprintf(stderr, "Wrote %lld %dd points to disk in %.4f seconds (%.2f points per millisecond)\n", n, d, t, ((n+0.0)/(t*1000)));

    t = -omp_get_wtime();
    PointStore B = PointStore::read_points("points.bin");
    t += omp_get_wtime();

    fprintf(stderr, "Read %lld %dd points from disk in %.4f seconds (%.2f points per millisecond)\n", n, d, t, ((n+0.0)/(t*1000)));
    fprintf(stderr, "read_points() %s\n", (A == B)? "succeeded" : "failed");


    return 0;
}
