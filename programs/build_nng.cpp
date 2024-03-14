#include <vector>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>
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
    double epsilon;
    char *infname = NULL;
    char *outfname = NULL;
    int nthreads = 1;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options] [..]\n", argv[0]);
        fprintf(stderr, "Options: -i FILE    input cover tree filename [required]\n");
        fprintf(stderr, "         -r FLOAT   epsilon radius [required]\n");
        fprintf(stderr, "         -o FILE    neighborhood graph file [default: none]\n");
        fprintf(stderr, "         -t INT     number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "         -v         log to stderr\n");
        fprintf(stderr, "         -h         help message\n");
        return 1;
    }

    infname = read_string_arg(argc, argv, "-i", NULL);
    epsilon = read_double_arg(argc, argv, "-r", NULL);
    outfname = read_string_arg(argc, argv, "-o", &outfname);
    nthreads = read_int_arg(argc, argv, "-t", &nthreads);
    verbose = !!(find_arg_idx(argc, argv, "-v") >= 0);

    omp_set_num_threads(std::min(nthreads, omp_get_max_threads()));
    fprintf(stderr, "Parameters: [radius=%.2f, infname='%s', outfname='%s', nthreads=%d]\n", epsilon, infname, outfname, nthreads);

    t = -omp_get_wtime();
    CoverTree tree;
    tree.read_from_file(infname);
    t += omp_get_wtime();

    fprintf(stderr, "Read cover tree from file '%s' in %.3f seconds\n", infname, t);
    if (verbose) tree.print_info(stderr);

    index_t n = tree.num_points();
    const float *p = tree.getdata();
    int d = tree.getdim();

    std::vector<std::vector<index_t>> graph;
    graph.resize(n);

    t = -omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 100)
    for (index_t i = 0; i < n; ++i) graph[i] = tree.radii_query(&p[d*i], epsilon);

    index_t m = 0;

    for (index_t i = 0; i < n; ++i)
    {
        m += graph[i].size();
        std::sort(graph[i].begin(), graph[i].end());
    }

    t += omp_get_wtime();

    fprintf(stderr, "Constructed neighborhood graph with %lld vertices and %lld edges with average degree %.2f in %.4f seconds\n", n, m, (m+0.0)/(n+0.0), t);

    if (outfname)
    {
        t = -omp_get_wtime();

        FILE *f = fopen(outfname, "w");

        fprintf(f, "%lld\t%lld\t%lld\n", n, n, m);

        for (index_t i = 0; i < n; ++i)
            for (index_t j : graph[i])
                fprintf(f, "%lld\t%lld\n", i, j);

        fclose(f);
        t += omp_get_wtime();

        fprintf(stderr, "Wrote graph to file '%s' in %.3f seconds\n", outfname, t);
    }

    return 0;
}
