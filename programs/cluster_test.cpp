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

std::vector<std::vector<index_t>>
get_clusters(const CoverTree& tree, const float *p, size_t n, int d, double epsilon);

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
        fprintf(stderr, "         -o FILE    cluster filename [default: none]\n");
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

    t = -omp_get_wtime();
    CoverTree tree;
    tree.read_from_file(infname);
    t += omp_get_wtime();

    if (verbose) fprintf(stderr, "Read cover tree from file '%s' in %.3f seconds\n", infname, t);
    if (verbose) tree.print_info();

    auto clusters = get_clusters(tree, tree.getdata(), tree.num_points(), tree.getdim(), epsilon);

    if (outfname)
    {
        std::ofstream stream(outfname);

        for (const auto& cluster : clusters)
        {
            std::copy(cluster.cbegin(), cluster.cend()-1, std::ostream_iterator<index_t>(stream, " "));
            stream << cluster.back() << std::endl;
        }

        stream.close();
    }

    return 0;
}

std::vector<std::vector<index_t>> connected_components(const std::vector<std::vector<index_t>>& graph)
{
    index_t n = graph.size();
    std::vector<std::vector<index_t>> components;
    std::vector<index_t> stack, component;
    std::vector<bool> visited(n, false);
    index_t s = 0;

    for (;;)
    {
        auto it = std::find(visited.begin() + s, visited.end(), false);

        if (it == visited.end()) break;
        s = it - visited.begin();

        component.clear();
        stack.resize(1);
        stack[0] = s;

        while (stack.size() != 0)
        {
            index_t u = stack.back(); stack.pop_back();

            if (visited[u])
                continue;

            visited[u] = true;
            component.push_back(u);

            std::copy_if(graph[u].begin(), graph[u].end(), std::back_inserter(stack), [&](index_t v) { return !visited[v]; });
        }

        std::sort(component.begin(), component.end());
        components.push_back(component);
    }

    return components;
}

std::vector<std::vector<index_t>>
get_clusters(const CoverTree& tree, const float *p, size_t n, int d, double epsilon)
{
    double t;

    std::vector<std::vector<index_t>> graph;
    graph.resize(n);

    t = -omp_get_wtime();
    #pragma omp parallel for schedule(dynamic, 100)
    for (index_t i = 0; i < n; ++i) graph[i] = tree.radii_query(&p[d*i], epsilon);
    t += omp_get_wtime();

    if (verbose) fprintf(stderr, "Constructed neighborhood graph in %.3f seconds\n", t);

    t = -omp_get_wtime();
    const auto& clusters = connected_components(graph);
    t += omp_get_wtime();

    if (verbose) fprintf(stderr, "Computed connected components in %.3f seconds\n", t);

    return clusters;
}
