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
#include "read_args.h"

static int verbose = 0;

namespace fs = std::filesystem;

std::vector<float> read_vecs_file(const char *fname, int *dim, size_t *npts);
void write_vecs_file(const char *fname, int dim, const std::vector<float>& pdata);

std::vector<std::vector<index_t>>
get_clusters(const CoverTree& tree, const float *p, size_t n, int d, double base, double epsilon);

int main(int argc, char *argv[])
{
    double t;
    double base = 2., epsilon = 1.;
    char *infname = NULL;
    char *outfname = NULL;
    int nthreads = 1;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options] [..]\n", argv[0]);
        fprintf(stderr, "Options: -i FILE    input filename [required]\n");
        fprintf(stderr, "         -r FLOAT   epsilon radius [default: %.2f]\n", epsilon);
        fprintf(stderr, "         -b FLOAT   cover tree base [default: %.2f]\n", base);
        fprintf(stderr, "         -o FILE    cluster filename [default: none]\n");
        fprintf(stderr, "         -t INT     number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "         -v         log to stderr\n");
        fprintf(stderr, "         -B         skip clustering step\n");
        fprintf(stderr, "         -h         help message\n");
        return 1;
    }

    base = read_double_arg(argc, argv, "-b", &base);
    epsilon = read_double_arg(argc, argv, "-r", &epsilon);
    outfname = read_string_arg(argc, argv, "-o", &outfname);
    infname = read_string_arg(argc, argv, "-i", NULL);
    nthreads = read_int_arg(argc, argv, "-t", &nthreads);
    verbose = !!(find_arg_idx(argc, argv, "-v") >= 0);

    omp_set_num_threads(std::min(nthreads, omp_get_max_threads()));

    int d;
    size_t n;

    t = -omp_get_wtime();
    auto pointmem = read_vecs_file(infname, &d, &n);
    t += omp_get_wtime();

    const float *p = pointmem.data();

    if (verbose) fprintf(stderr, "Read points file '%s' (%lld points of dimension %d) in %.3f seconds\n", infname, n, d, t);

    t = -omp_get_wtime();
    CoverTree tree(p, n, d, base);
    t += omp_get_wtime();

    if (verbose) fprintf(stderr, "Constructed cover tree in %.3f seconds\n", t);

    assert(tree.is_full());
    assert(tree.is_nested());
    assert(tree.is_covering());

    if (verbose) tree.print_info();

    if (find_arg_idx(argc, argv, "-B") >= 0)
        return 0;

    auto clusters = get_clusters(tree, p, n, d, base, epsilon);

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

void write_vecs_file(const char *fname, int dim, const std::vector<float>& pdata)
{
    size_t n = pdata.size() / dim;
    const float *p = pdata.data();
    FILE *f;

    f = fopen(fname, "wb");

    for (size_t i = 0; i < n; ++i)
    {
        fwrite(&dim, sizeof(int), 1, f);
        fwrite(&p[i*dim], sizeof(float), dim, f);
    }

    fclose(f);
}

std::vector<float> read_vecs_file(const char *fname, int *dim, size_t *npts)
{
    std::vector<float> pdata;
    size_t filesize, n;
    int d;
    FILE *f;
    fs::path path = fname;

    if (!fs::exists(path) || !fs::is_regular_file(path))
    {
        std::cerr << "error: unable to open " << std::quoted(fname) << std::endl;
        exit(1);
    }

    filesize = fs::file_size(path);

    f = fopen(fname, "rb");
    fread(&d, sizeof(int), 1, f);
    n = filesize / (4*(d+1));
    fseek(f, SEEK_SET, 0);

    if (n == 0 || d <= 0)
    {
        std::cerr << "error: unable to open " << std::quoted(fname) << std::endl;
        exit(1);
    }

    pdata.resize(n*d);
    float *p = pdata.data();

    for (size_t i = 0; i < n; ++i)
    {
        fread(&d, sizeof(int), 1, f);
        fread(&p[i*d], sizeof(float), d, f);
    }

    fclose(f);

    if (dim) *dim = d;
    if (npts) *npts = n;

    return pdata;
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
get_clusters(const CoverTree& tree, const float *p, size_t n, int d, double base, double epsilon)
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
