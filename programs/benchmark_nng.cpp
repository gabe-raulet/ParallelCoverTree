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
#include "CoverTree.h"
#include "read_args.h"

static int max_threads;

extern double distance(const float *p, const float *q, int d);

void run_test(size_t n, int d, double var, double radius);
std::vector<std::vector<int64_t>> build_nng(const CoverTree& tree, double radius, size_t *n_edges);
std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, size_t *n_edges);

int main(int argc, char *argv[])
{
    max_threads = omp_get_max_threads();

    int64_t npoints;
    int dim;
    double radius;
    double var;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [..]\n", argv[0]);
        fprintf(stderr, "Parameters: -n INT    number of points\n");
        fprintf(stderr, "            -d INT    dimension\n");
        fprintf(stderr, "            -V FLOAT  variance\n");
        fprintf(stderr, "            -r FLOAT  epsilon radius\n");
        fprintf(stderr, "            -h        help message\n");
        return 1;
    }

    npoints = read_formatted_int_arg(argc, argv, "-n", NULL);
    dim = read_int_arg(argc, argv, "-d", NULL);
    radius = read_double_arg(argc, argv, "-r", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);

    run_test(static_cast<size_t>(npoints), dim, var, radius);
}

std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, size_t *n_edges)
{
    size_t n = pointmem.size() / dim;
    const float *pdata = pointmem.data();
    std::vector<std::vector<int64_t>> graph(n);

    size_t m = 0;

    for (size_t u = 0; u < n; ++u)
    {
        graph[u].push_back(u);

        for (size_t v = u+1; v < n; ++v)
            if (distance(&pdata[dim*u], &pdata[dim*v], dim) <= radius)
            {
                graph[u].push_back(v);
                graph[v].push_back(u);
            }

        m += graph[u].size();
    }

    *n_edges = m;
    return graph;
}

std::vector<std::vector<int64_t>> build_nng(const CoverTree& tree, double radius, size_t *n_edges)
{
    size_t m = 0;
    size_t n = tree.num_points();
    const float *p = tree.getdata();
    int d = tree.getdim();

    std::vector<std::vector<int64_t>> graph(n);

    #pragma omp parallel for schedule(dynamic, 100) reduction(+:m)
    for (index_t i = 0; i < n; ++i) 
    {
        graph[i] = tree.radii_query(&p[d*i], radius);
        m += graph[i].size();
    }
    
    *n_edges = m;
    return graph;
}

void run_test(size_t n, int d, double var, double radius)
{
    double t;

    t = -omp_get_wtime();
    std::random_device rd;
    std::default_random_engine gen;
    std::normal_distribution dis{0.0, std::sqrt(var)};
    std::vector<float> pointmem(d*n);
    std::generate(pointmem.begin(), pointmem.end(), [&]() { return dis(gen); });
    t += omp_get_wtime();

    fprintf(stderr, "(generate points) :: [n=%lld,d=%d,var=%.2f,radius=%.2f] :: [%.4f (seconds)]\n", n, d, var, radius, t);

    const float *p = pointmem.data();
    double base = 2.;

    t = -omp_get_wtime();
    CoverTree tree(p, n, d, base);
    t += omp_get_wtime();

    fprintf(stderr, "(build cover tree) :: [n_verts=%lld,base=%.2f] :: [%.4f (seconds)]\n", tree.num_vertices(), base, t);

    size_t m;

    t = -omp_get_wtime();
    const auto& graph = build_nng_bf(pointmem, d, radius, &m);
    t += omp_get_wtime();

    fprintf(stderr, "(build nng bf) :: [n_verts=%lld,n_edges=%lld,avg_deg=%.2f] :: [%.4f (seconds)]\n", n, m, (m+0.0)/(n+0.0), t);

    for (int nthreads = 1; nthreads <= max_threads; nthreads *= 2)
    {
        omp_set_num_threads(nthreads);

        t = -omp_get_wtime();
        const auto& graph = build_nng(tree, radius, &m);
        t += omp_get_wtime();

        fprintf(stderr, "(build nng) :: [n_verts=%lld,n_edges=%lld,avg_deg=%.2f] :: [nthreads=%d] :: [%.4f (seconds)]\n", n, m, (m+0.0)/(n+0.0), nthreads, t);
    }
}
