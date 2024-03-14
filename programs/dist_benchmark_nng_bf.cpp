#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <random>
#include <limits.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <mpi.h>
#include "read_args.h"

double distance(const float *p, const float *q, int d);
std::vector<float> generate_points(int64_t n, int d, double var, int seed, MPI_Comm comm);
std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, int64_t *n_edges, MPI_Comm comm);
void write_degrees_to_file(const char *fname, std::vector<std::vector<int64_t>>& graph, MPI_Comm comm);

int main(int argc, char *argv[])
{
    int64_t npoints;
    int dim;
    double radius;
    double var;
    int seed = -1;
    char *outfname = NULL;

    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        if (myrank == 0)
        {
            fprintf(stderr, "Usage: %s [..]\n", argv[0]);
            fprintf(stderr, "Parameters: -n INT    number of points\n");
            fprintf(stderr, "            -d INT    dimension\n");
            fprintf(stderr, "            -V FLOAT  variance\n");
            fprintf(stderr, "            -r FLOAT  epsilon radius\n");
            fprintf(stderr, "            -o FILE   graph degrees file [optional]\n");
            fprintf(stderr, "            -s INT    rng seed [default: random]\n");
            fprintf(stderr, "            -h        help message\n");
        }

        MPI_Finalize();
        return 0;
    }

    npoints = read_formatted_int_arg(argc, argv, "-n", NULL);
    dim = read_int_arg(argc, argv, "-d", NULL);
    radius = read_double_arg(argc, argv, "-r", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);
    seed = read_int_arg(argc, argv, "-s", &seed);

    if (find_arg_idx(argc, argv, "-o") >= 0)
    {
        outfname = read_string_arg(argc, argv, "-o", NULL);
    }

    double t, maxtime, proctime;

    MPI_Barrier(MPI_COMM_WORLD);
    t = -MPI_Wtime();

    auto pointmem = generate_points(npoints, dim, var, seed, MPI_COMM_WORLD);

    t += MPI_Wtime();
    MPI_Reduce(&t, &maxtime,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myrank == 0) fprintf(stderr, "(generate_points) :: [n_points=%lld,dim=%d,var=%.2f] :: [maxtime=%.4f,avgtime=%.4f (seconds)]\n", npoints, dim, var, maxtime, proctime/nprocs);

    int64_t m, n_edges;

    MPI_Barrier(MPI_COMM_WORLD);
    t = -MPI_Wtime();

    auto local_dist_graph = build_nng_bf(pointmem, dim, radius, &m, MPI_COMM_WORLD);

    MPI_Reduce(&m, &n_edges, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    t += MPI_Wtime();
    MPI_Reduce(&t, &maxtime,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myrank == 0) fprintf(stderr, "(build_nng_bf) :: [n_verts=%lld,n_edges=%lld,avg_deg=%.2f,radius=%.2f] :: [maxtime=%.4f,avgtime=%.4f (seconds)]\n", npoints, n_edges, (n_edges+0.0)/npoints, radius, maxtime, proctime/nprocs);

    if (outfname)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        t = -MPI_Wtime();

        write_degrees_to_file(outfname, local_dist_graph, MPI_COMM_WORLD);

        t += MPI_Wtime();
        MPI_Reduce(&t, &maxtime,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (myrank == 0) fprintf(stderr, "(write_degrees_to_file) :: [outfname='%s'] :: [maxtime=%.4f,avgtime=%.4f (seconds)]\n", outfname, maxtime, proctime/nprocs);
    }

    MPI_Finalize();
    return 0;
}

std::vector<float> generate_points(int64_t n, int d, double var, int seed, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    std::vector<float> pointmem(d*n);

    if (myrank == 0)
    {
        std::random_device rd;
        std::default_random_engine gen(seed < 0? rd() : seed*17);
        std::normal_distribution dis{0.0, std::sqrt(var)};
        std::generate(pointmem.begin(), pointmem.end(), [&]() { return dis(gen); });
    }

    MPI_Bcast(pointmem.data(), static_cast<int>(n*d), MPI_FLOAT, 0, comm);
    return std::move(pointmem);
}

std::vector<std::vector<int64_t>> build_nng_bf(const std::vector<float>& pointmem, int dim, double radius, int64_t *n_edges, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int64_t n = pointmem.size() / dim;
    const float *pdata = pointmem.data();

    int64_t chunksize = n / nprocs;
    int64_t start = myrank * chunksize;
    int64_t stop = std::min(start + chunksize, n);
    int64_t mychunksize = stop - start;

    std::vector<std::vector<int64_t>> graph(mychunksize);

    int64_t m = 0;

    for (int64_t u = start; u < stop; ++u)
    {
        for (int64_t v = 0; v < n; ++v)
            if (distance(&pdata[dim*u], &pdata[dim*v], dim) <= radius)
                graph[u-start].push_back(v);

        m += graph[u-start].size();
    }

    *n_edges = m;

    return std::move(graph);
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

void write_degrees_to_file(const char *fname, std::vector<std::vector<int64_t>>& graph, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int64_t my_n = graph.size();
    std::vector<int64_t> mydegrees(my_n);

    for (int64_t i = 0; i < my_n; ++i)
        mydegrees[i] = graph[i].size();

    std::vector<int> recvcounts, displs;

    if (myrank == 0) recvcounts.resize(nprocs);

    int sendcount = my_n;
    MPI_Gather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, comm);

    std::vector<int64_t> degrees;

    if (myrank == 0)
    {
        displs.resize(nprocs);
        displs[0] = 0;
        std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
        degrees.resize(displs.back() + recvcounts.back());
    }

    MPI_Gatherv(mydegrees.data(), sendcount, MPI_INT64_T, degrees.data(), recvcounts.data(), displs.data(), MPI_INT64_T, 0, comm);

    if (myrank == 0)
    {
        FILE *f = fopen(fname, "w");
        for (auto d : degrees) fprintf(f, "%lld\n", d);
        fclose(f);
    }
}
