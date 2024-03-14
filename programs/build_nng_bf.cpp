#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "VectorIO.h"
#include "read_args.h"

double distance(const float *p, const float *q, int d);

int main(int argc, char *argv[])
{
    double radius, t;
    char *infname, *outfname;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -r FLOAT  epsilon radius [required]\n");
        fprintf(stderr, "         -i FILE   input filename [required]\n");
        fprintf(stderr, "         -o FILE   output filename [required]\n");
        return 1;
    }

    radius = read_double_arg(argc, argv, "-r", NULL);
    infname = read_string_arg(argc, argv, "-i", NULL);
    outfname = read_string_arg(argc, argv, "-o", NULL);

    fprintf(stderr, "Parameters: [radius=%.2f, infname='%s', outfname='%s']\n", radius, infname, outfname);

    int d;
    size_t n;
    std::vector<float> pointmem;

    t = -omp_get_wtime();
    pointmem = read_vecs_file(infname, &d, &n);
    t += omp_get_wtime();

    fprintf(stderr, "Read %lld points of dimension %d from file '%s' in %.4f seconds\n", n, d, infname, t);

    t = -omp_get_wtime();

    const float *pdata = pointmem.data();
    std::vector<std::vector<size_t>> graph(n);

    size_t m = 0;

    for (size_t u = 0; u < n; ++u)
    {
        graph[u].push_back(u);

        for (size_t v = u+1; v < n; ++v)
            if (distance(&pdata[d*u], &pdata[d*v], d) <= radius)
            {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }

        m += graph[u].size();
    }

    t += omp_get_wtime();

    fprintf(stderr, "Constructed neighborhood graph with %lld vertices and %lld edges with average degree %.2f in %.4f seconds\n", n, m, (m+0.0)/(n+0.0), t);

    t = -omp_get_wtime();

    FILE *f = fopen(outfname, "w");

    fprintf(f, "%lld\t%lld\t%lld\n", n, n, m);

    for (size_t u = 0; u < n; ++u)
        for (size_t v : graph[u])
            fprintf(f, "%lld\t%lld\n", u, v);

    fclose(f);

    t += omp_get_wtime();

    fprintf(stderr, "Wrote graph to file '%s' in %.4f seconds\n", outfname, t);

    return 0;
}

double distance(const float *p, const float *q, int d)
{
    double sum = 0., di;

    for (int i = 0; i < d; ++i)
    {
        di = p[i] - q[i];
        sum += di*di;
    }

    return std::sqrt(sum);
}
