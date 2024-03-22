#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdio.h>
#include <mpi.h>
#include "MPITimer.h"
#include "read_args.h"

double distance(const float *p, const float *q, int d);
std::vector<float> generate_points(int64_t n, int d, double var, int seed);
std::vector<float> scatter_points(const std::vector<float>& points, int d, int root, MPI_Comm comm);
std::vector<float> gather_points(const std::vector<float>& mypoints, int d, int root, MPI_Comm comm);
std::vector<int64_t> serial_greedy_permutation(const std::vector<float>& points, int d, int64_t start);
std::vector<int64_t> parallel_greedy_permutation(const std::vector<float>& mypoints, int d, int64_t start, int root, MPI_Comm comm);
void serial_greedy_update(std::vector<double>& dists, const std::vector<float>& points, int d, const float *cur);
void parallel_greedy_update(std::vector<double>& mydists, const std::vector<float>& mypoints, int d, const float *curpoint, int currank, MPI_Comm comm);
int64_t serial_argmax(const std::vector<double>& v);
int64_t parallel_argmax(const std::vector<double>& v, MPI_Comm comm);

int main(int argc, char *argv[])
{
    int64_t n;
    int d;
    double var;
    int seed = -1;

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
            fprintf(stderr, "            -s INT    rng seed [default: random]\n");
            fprintf(stderr, "            -h        help message\n");
        }

        MPI_Finalize();
        return 0;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    d = read_int_arg(argc, argv, "-d", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);
    seed = read_int_arg(argc, argv, "-s", &seed);

    assert(n >= 1);
    MPITimer timer(MPI_COMM_WORLD, 0);

    double serial_time, parallel_time;
    std::vector<float> points, mypoints;
    std::vector<int64_t> perm, perm2;

    timer.start_timer();
    if (!myrank) points = generate_points(n, d, var, seed);
    timer.stop_timer();

    if (!myrank) fprintf(stderr, "(generate_points) :: %lld points :: [maxtime=%.4f, avgtime=%.4f (seconds)]\n", n, timer.get_max_time(), timer.get_avg_time());

    timer.start_timer();
    mypoints = scatter_points(points, d, 0, MPI_COMM_WORLD);
    timer.stop_timer();

    if (!myrank) fprintf(stderr, "(scatter_points) :: %lld points :: [maxtime=%.4f, avgtime=%.4f (seconds)]\n", n, timer.get_max_time(), timer.get_avg_time());

    timer.start_timer();
    if (!myrank) perm = serial_greedy_permutation(points, d, 0);
    timer.stop_timer();
    serial_time = timer.get_max_time();

    //if (!myrank) std::copy(perm.begin(), perm.end(), std::ostream_iterator<int64_t>(std::cout, "\n"));

    if (!myrank) fprintf(stderr, "(serial_greedy_permutation) :: [maxtime=%.4f, avgtime=%.4f (seconds)]\n", timer.get_max_time(), timer.get_avg_time());

    timer.start_timer();
    perm2 = parallel_greedy_permutation(mypoints, d, 0, 0, MPI_COMM_WORLD);
    timer.stop_timer();
    parallel_time = timer.get_max_time();

    if (!myrank) fprintf(stderr, "(parallel_greedy_permutation) :: [maxtime=%.4f, avgtime=%.4f (seconds)]\n", timer.get_max_time(), timer.get_avg_time());

    if (!myrank)
    {
        if (perm == perm2) fprintf(stderr, "Successfully computed greedy permutation in parallel\n");
        else fprintf(stderr, "Failed trying to compute greedy permutation in parallel\n");

        fprintf(stderr, "%.4f speedup for parallel algorithm\n", serial_time / parallel_time);
    }


    MPI_Finalize();
    return 0;
}

std::vector<float> generate_points(int64_t n, int d, double var, int seed)
{
    std::vector<float> pointmem(d*n);
    std::random_device rd;
    std::default_random_engine gen(seed < 0? rd() : seed*17);
    std::normal_distribution dis{0.0, std::sqrt(var)};
    std::generate(pointmem.begin(), pointmem.end(), [&]() { return dis(gen); });
    return std::move(pointmem);
}

std::vector<float> scatter_points(const std::vector<float>& points, int d, int root, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    std::vector<int> sendcounts, displs;
    std::vector<float> mypoints;

    if (myrank == root)
    {
        sendcounts.resize(nprocs);
        displs.resize(nprocs);

        int n = points.size() / d;
        std::fill(sendcounts.begin(), sendcounts.end(), n/nprocs);
        sendcounts.back() = n - (nprocs-1)*(n/nprocs);
        std::for_each(sendcounts.begin(), sendcounts.end(), [&](auto& item) { item *= d; });

        displs[0] = 0;
        std::partial_sum(sendcounts.begin(), sendcounts.end()-1, displs.begin()+1);
    }

    int mysize;
    MPI_Scatter(sendcounts.data(), 1, MPI_INT, &mysize, 1, MPI_INT, root, comm);

    mypoints.resize(mysize);
    MPI_Scatterv(points.data(), sendcounts.data(), displs.data(), MPI_FLOAT, mypoints.data(), mysize, MPI_FLOAT, root, comm);

    return mypoints;
}

std::vector<float> gather_points(const std::vector<float>& mypoints, int d, int root, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    std::vector<float> points;
    std::vector<int> recvcounts, displs;

    if (myrank == root)
    {
        recvcounts.resize(nprocs);
        displs.resize(nprocs);
    }

    int mysize = mypoints.size();
    MPI_Gather(&mysize, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, root, comm);

    if (myrank == root)
    {
        displs[0] = 0;
        std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
        points.resize(displs.back() + recvcounts.back());
    }

    MPI_Gatherv(mypoints.data(), mysize, MPI_FLOAT, points.data(), recvcounts.data(), displs.data(), MPI_FLOAT, root, comm);
    return points;
}

void serial_greedy_update(std::vector<double>& dists, const std::vector<float>& points, int d, const float *cur)
{
    int64_t n = dists.size();
    assert(n == points.size() / d);
    const float *p = &points[0];

    for (int64_t j = 0; j < n; ++j, p += d)
    {
        dists[j] = std::min(dists[j], distance(p, cur, d));
    }
}

int64_t serial_argmax(const std::vector<double>& v)
{
    return std::distance(v.begin(), std::max_element(v.begin(), v.end()));
}

std::vector<int64_t> serial_greedy_permutation(const std::vector<float>& points, int d, int64_t start)
{
    int64_t n = points.size() / d;
    std::vector<int64_t> perm;
    std::vector<double> dists(n);
    perm.reserve(n);
    const float *curpoint;

    int64_t cur = start;
    std::fill(dists.begin(), dists.end(), std::numeric_limits<double>::max());

    for (int64_t i = 0; i < n; ++i)
    {
        perm.push_back(cur);
        curpoint = &points[cur*d];
        serial_greedy_update(dists, points, d, curpoint); /* updates dists */
        cur = serial_argmax(dists);
    }

    return perm;
}

int get_cur_rank(int64_t mysize, int64_t index, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int64_t myoffset;
    MPI_Exscan(&mysize, &myoffset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (myrank == 0) myoffset = 0;

    std::vector<int> owners(nprocs, 0);

    if (myoffset <= index && index < myoffset + mysize)
        owners[myrank] = 1;

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, owners.data(), 1, MPI_INT, comm);

    return std::find(owners.begin(), owners.end(), 1) - owners.begin();
}

void parallel_greedy_update(std::vector<double>& mydists, const std::vector<float>& mypoints, int d, const float *curpoint, int currank, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    std::vector<float> cur(d);
    if (myrank == currank) std::copy(curpoint, curpoint + d, cur.begin());
    MPI_Bcast(cur.data(), d, MPI_FLOAT, currank, comm);

    int64_t mysize = mydists.size();
    const float *p = &mypoints[0];

    for (int64_t j = 0; j < mysize; ++j, p += d)
    {
        mydists[j] = std::min(mydists[j], distance(p, cur.data(), d));
    }
}

int64_t parallel_argmax(const std::vector<double>& v, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int64_t mysize = v.size();
    int64_t myoffset;
    MPI_Exscan(&mysize, &myoffset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (myrank == 0) myoffset = 0;

    int64_t myargmax = serial_argmax(v);
    double mymax = v[myargmax];

    std::vector<double> maxes(nprocs);
    std::vector<int64_t> argmaxes(nprocs);
    maxes[myrank] = mymax;
    argmaxes[myrank] = myargmax + myoffset;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_DOUBLE, maxes.data(), 1, MPI_DOUBLE, comm);
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT64_T, argmaxes.data(), 1, MPI_INT64_T, comm);

    int maxrank = serial_argmax(maxes);
    return argmaxes[maxrank];
}

std::vector<int64_t> parallel_greedy_permutation(const std::vector<float>& mypoints, int d, int64_t start, int root, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int64_t mysize = mypoints.size() / d;
    int64_t myoffset, n;

    MPI_Exscan(&mysize, &myoffset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (myrank == 0) myoffset = 0;

    MPI_Allreduce(&mysize, &n, 1, MPI_INT64_T, MPI_SUM, comm);

    std::vector<int64_t> perm;
    std::vector<double> mydists(mysize);
    const float *curpoint;

    int64_t cur = start;
    std::fill(mydists.begin(), mydists.end(), std::numeric_limits<double>::max());

    for (int64_t i = 0; i < n; ++i)
    {
        if (myrank == root) perm.push_back(cur);
        int currank = get_cur_rank(mydists.size(), cur, comm);
        curpoint = currank == myrank? &mypoints[(cur-myoffset)*d] : NULL;
        parallel_greedy_update(mydists, mypoints, d, curpoint, currank, comm);
        cur = parallel_argmax(mydists, comm);
    }

    return perm;
}

double distance(const float *p, const float *q, int d)
{
    double sum = 0.0, delta;

    for (int i = 0; i < d; ++i)
    {
        delta = p[i] - q[i];
        sum += delta * delta;
    }

    return std::sqrt(sum);
}
