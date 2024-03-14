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
static int nthreads = 1;

double distance(const float *p, const float *q, int d);
std::vector<float> generate_points(size_t n, int d, double var, int seed);
std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, size_t *n_edges);
void write_degrees_to_file(const char *fname, std::vector<std::vector<int64_t>>& graph);

int main(int argc, char *argv[])
{
    max_threads = omp_get_max_threads();

    int64_t npoints;
    int dim;
    double radius;
    double var;
    int seed = -1;
    char *outfname = NULL;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [..]\n", argv[0]);
        fprintf(stderr, "Parameters: -n INT    number of points\n");
        fprintf(stderr, "            -d INT    dimension\n");
        fprintf(stderr, "            -V FLOAT  variance\n");
        fprintf(stderr, "            -r FLOAT  epsilon radius\n");
        fprintf(stderr, "            -o FILE   graph degrees file [optional]\n");
        fprintf(stderr, "            -t INT    number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "            -s INT    rng seed [default: random]\n");
        fprintf(stderr, "            -h        help message\n");
        return 1;
    }

    npoints = read_formatted_int_arg(argc, argv, "-n", NULL);
    dim = read_int_arg(argc, argv, "-d", NULL);
    radius = read_double_arg(argc, argv, "-r", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);
    seed = read_int_arg(argc, argv, "-s", &seed);
    nthreads = std::min(read_int_arg(argc, argv, "-t", &nthreads), max_threads);

    if (find_arg_idx(argc, argv, "-o") >= 0)
    {
        outfname = read_string_arg(argc, argv, "-o", NULL);
    }

    double t;

    t = -omp_get_wtime();
    auto pointmem = generate_points(static_cast<size_t>(npoints), dim, var, seed);
    t += omp_get_wtime();

    fprintf(stderr, "(generate_points) :: [n_points=%lld,dim=%d,var=%.2f] :: [%.4f seconds]\n", npoints, dim, var, t);

    size_t m;

    t = -omp_get_wtime();
    auto graph = build_nng_bf(pointmem, dim, radius, &m);
    t += omp_get_wtime();

    fprintf(stderr, "(build_nng_bf) :: [n_verts=%lld,n_edges=%lld,avg_deg=%.2f,radius=%.2f,nthreads=%d] :: [%.4f seconds]\n", npoints, m, (m+0.0)/npoints, radius, nthreads, t);

    if (outfname)
    {
        t = -omp_get_wtime();
        write_degrees_to_file(outfname, graph);
        t += omp_get_wtime();

        fprintf(stderr, "(write_degrees_to_file) :: [outfname='%s',nthreads=%d] :: [%.4f seconds]\n", outfname, nthreads, t);
    }

    return 0;
}

std::vector<float> generate_points(size_t n, int d, double var, int seed)
{
    std::random_device rd;
    std::default_random_engine gen(seed < 0? rd() : seed*17);
    std::normal_distribution dis{0.0, std::sqrt(var)};
    std::vector<float> pointmem(d*n);
    std::generate(pointmem.begin(), pointmem.end(), [&]() { return dis(gen); });
    return std::move(pointmem);
}

std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, size_t *n_edges)
{
    omp_set_num_threads(nthreads);

    size_t n = pointmem.size() / dim;
    const float *pdata = pointmem.data();
    std::vector<std::vector<int64_t>> graph(n);

    size_t m = 0;

    #pragma omp parallel for reduction(+:m)
    for (size_t u = 0; u < n; ++u)
    {
        for (size_t v = 0; v < n; ++v)
            if (distance(&pdata[dim*u], &pdata[dim*v], dim) <= radius)
                graph[u].push_back(static_cast<int64_t>(v));

        m += graph[u].size();
    }

    *n_edges = m;
    return std::move(graph);
}

void write_degrees_to_file(const char *fname, std::vector<std::vector<int64_t>>& graph)
{
    omp_set_num_threads(nthreads);

    size_t m = 0;
    size_t n = graph.size();

    std::vector<int64_t> degrees(n);

    #pragma omp parallel for reduction(+:m)
    for (size_t i = 0; i < n; ++i)
    {
        degrees[i] = graph[i].size();
        m += degrees[i];
    }

    FILE *f = fopen(fname, "w");

    for (size_t i = 0; i < n; ++i)
    {
        fprintf(f, "%lld\n", degrees[i]);
    }

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
