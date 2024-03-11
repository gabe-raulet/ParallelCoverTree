#include <vector>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <sys/stat.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "CoverTree.h"
#include "VectorIO.h"
#include "read_args.h"

static int verbose = 0;

int main(int argc, char *argv[])
{
    double t;
    double base = 2.;
    char *infname = NULL, *outfname = NULL;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options] [..]\n", argv[0]);
        fprintf(stderr, "Options: -i FILE    input vector filename [required]\n");
        fprintf(stderr, "         -o FILE    output cover tree filename [required]\n");
        fprintf(stderr, "         -b FLOAT   cover tree base [default: %.2f]\n", base);
        fprintf(stderr, "         -v         log to stderr\n");
        fprintf(stderr, "         -h         help message\n");
        return 1;
    }

    infname = read_string_arg(argc, argv, "-i", NULL);
    outfname = read_string_arg(argc, argv, "-o", NULL);
    base = read_double_arg(argc, argv, "-b", &base);
    verbose = !!(find_arg_idx(argc, argv, "-v") >= 0);

    int d;
    size_t n;

    t = -omp_get_wtime();
    auto pointmem = read_vecs_file(infname, &d, &n);
    t += omp_get_wtime();

    fprintf(stderr, "Read points file '%s' (%lld points of dimension %d) in %.3f seconds\n", infname, n, d, t);

    const float *p = pointmem.data();

    t = -omp_get_wtime();
    CoverTree tree(p, n, d, base);
    t += omp_get_wtime();

    fprintf(stderr, "Constructed cover tree in %.3f seconds\n", t);

    assert(tree.is_full());
    assert(tree.is_nested());
    assert(tree.is_covering());

    if (verbose) tree.print_info();

    t = -omp_get_wtime();
    tree.write_to_file(outfname);
    t += omp_get_wtime();

    fprintf(stderr, "Wrote cover tree to file '%s' in %.3f seconds\n", outfname, t);

    return 0;
}
