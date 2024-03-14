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
#include "read_args.h"

static int max_threads;

double distance(const float *p, const float *q, int d);
std::vector<float> generate_points(size_t n, int d, double var);
std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, size_t *n_edges);
void write_graph_to_file(const char *outfname, std::vector<std::vector<int64_t>>& graph);

int main(int argc, char *argv[])
{
    max_threads = omp_get_max_threads();

    int64_t npoints;
    int dim;
    double radius;
    double var;
    int nthreads = 1;
    char *outfname = NULL;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [..]\n", argv[0]);
        fprintf(stderr, "Parameters: -n INT    number of points\n");
        fprintf(stderr, "            -d INT    dimension\n");
        fprintf(stderr, "            -V FLOAT  variance\n");
        fprintf(stderr, "            -r FLOAT  epsilon radius\n");
        fprintf(stderr, "            -o FILE   output graph file [optional]\n");
        fprintf(stderr, "            -t INT    number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "            -h        help message\n");
        return 1;
    }

    npoints = read_formatted_int_arg(argc, argv, "-n", NULL);
    dim = read_int_arg(argc, argv, "-d", NULL);
    radius = read_double_arg(argc, argv, "-r", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);
    nthreads = std::min(read_int_arg(argc, argv, "-t", &nthreads), max_threads);

    if (find_arg_idx(argc, argv, "-o") >= 0)
    {
        outfname = read_string_arg(argc, argv, "-o", NULL);
    }

    double t;

    t = -omp_get_wtime();
    auto pointmem = generate_points(static_cast<size_t>(npoints), dim, var);
    t += omp_get_wtime();

    fprintf(stderr, "(generate_points) :: [n=%lld,d=%d,var=%.2f] :: [%.4f seconds]\n", npoints, dim, var, t);

    size_t m;

    t = -omp_get_wtime();
    auto graph = build_nng_bf(pointmem, dim, radius, &m);
    t += omp_get_wtime();

    fprintf(stderr, "(build_nng_bf) :: [n_edges=%lld,avg_deg=%.2f,radius=%.2f,nthreads=%d] :: [%.4f seconds]\n", m, (m+0.0)/npoints, radius, nthreads, t);

    if (outfname)
    {
        t = -omp_get_wtime();
        write_graph_to_file(outfname, graph);
        t += omp_get_wtime();

        fprintf(stderr, "(write_graph_to_file) :: [outfname='%s'] :: [%.4f seconds]\n", outfname, t);
    }

    return 0;
}

std::vector<float> generate_points(size_t n, int d, double var)
{
    std::random_device rd;
    std::default_random_engine gen;
    std::normal_distribution dis{0.0, std::sqrt(var)};
    std::vector<float> pointmem(d*n);
    std::generate(pointmem.begin(), pointmem.end(), [&]() { return dis(gen); });
    return std::move(pointmem);
}

std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, size_t *n_edges)
{
    size_t n = pointmem.size() / dim;
    const float *pdata = pointmem.data();
    std::vector<std::vector<int64_t>> graph(n);

    size_t m = 0;

    for (size_t u = 0; u < n; ++u)
    {
        for (size_t v = u+1; v < n; ++v)
            if (distance(&pdata[dim*u], &pdata[dim*v], dim) <= radius)
                graph[u].push_back(static_cast<int64_t>(v));

        m += graph[u].size();
    }

    *n_edges = m;
    return std::move(graph);
}

void write_graph_to_file(const char *outfname, std::vector<std::vector<int64_t>>& graph)
{
    size_t m = 0;
    size_t n = graph.size();

    for (auto& adj : graph)
    {
        m += adj.size();
        std::sort(adj.begin(), adj.end());
    }

    FILE *f = fopen(outfname, "w");

    fprintf(f, "%lld\t%lld\t%lld\n", n, n, m);

    for (int64_t u = 0; u < n; ++u)
        for (int64_t v : graph[u])
            fprintf(f, "%lld\t%lld\n", u, v);

    fclose(f);
}

double distance(const float *p, const float *q, int d)
{
    double sum = 0.0, di;

    for (int i = 0; i < d; ++i)
    {
        di = p[i] - q[i];
        sum += di*di;
    }

    return std::sqrt(sum);
}
