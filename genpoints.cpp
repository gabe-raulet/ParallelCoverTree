#include <iostream>
#include <random>
#include <vector>
#include <unistd.h>
#include <mpi.h>
#include "mpi_env.h"
#include "genpoints.h"
#include "ptraits.h"
#include "misc.h"

using namespace std;

using Comm = MPIEnv::Comm;

using PointTraits = SelectPoint<FPSIZE, PTDIM>::Traits;
using Real = PointTraits::Real;
using Point = PointTraits::Point;

template <integral Index>
void read_options(int argc, char *argv[], Index& n, char *&fname, double& var, int& seed, int& blocks, Comm comm);

void generate_points(vector<Point>& mypoints, size_t totsize, double var, int seed, int blocks, int root, Comm comm);
void write_points_file(const vector<Point>& mypoints, const char *fname, Comm comm);

int main_mpi(int argc, char *argv[]);
int main(int argc, char *argv[])
{
    MPIEnv::initialize(&argc, &argv);
    int ret = main_mpi(argc, argv);
    MPIEnv::finalize();
    return ret;
}

int main_mpi(int argc, char *argv[])
{
    Comm comm = Comm::comm_world();
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();
    timer.start_timer();

    size_t n;
    char *fname;
    double var = 10.0;
    int seed = -1;
    int blocks = 12;
    vector<Point> mypoints;

    read_options(argc, argv, n, fname, var, seed, blocks, comm);
    generate_points(mypoints, n, var, seed, blocks, 0, comm);
    write_points_file(mypoints, fname, comm);

    timer.stop_timer();

    if (!myrank) main_msg(argc, argv, timer.get_max_time(), timer.get_sum_time(), nprocs);

    return 0;
}

template <integral Index>
void read_options(int argc, char *argv[], Index& n, char *&fname, double& var, int& seed, int& blocks, Comm comm)
{
    double t, maxtime, sumtime;
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();
    timer.start_timer();

    auto usage = [&argv, &var, &seed, &blocks] (int err, bool isroot)
    {
        if (isroot)
        {
            fprintf(stderr, "Usage: %s [options] <numpts> <outfname>\n", argv[0]);
            fprintf(stderr, "Options: -V FLOAT variance [%.2f]\n", var);
            fprintf(stderr, "         -s INT   random seed [%s]\n", seed < 0? "random" : to_string(seed).c_str());
            fprintf(stderr, "         -b INT   number of seed blocks [%d]\n", blocks);
            fprintf(stderr, "         -h       help message\n");
        }

        MPIEnv::exit(err);
    };

    int c;
    while ((c = getopt(argc, argv, "V:s:b:h")) >= 0)
    {
        if      (c == 'V') var = atof(optarg);
        else if (c == 's') seed = atoi(optarg);
        else if (c == 'b') blocks = atoi(optarg);
        else if (c == 'h') usage(0, !myrank);
    }

    if (argc - optind < 2)
    {
        if (!myrank) fprintf(stderr, "[err::%s] missing argument(s)\n", __func__);
        usage(1, !myrank);
    }

    if (blocks < nprocs)
    {
        if (!myrank) fprintf(stderr, "[err::%s] blocks(%d) must be >= than nprocs(%d)\n", __func__, blocks, nprocs);
        usage(1, !myrank);
    }

    n = read_integer<Index>(argv[optind++]);
    fname = argv[optind];

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] :: [n=%lu,fname='%s',var=%.2f,seed=%d]\n", timer.get_max_time(), timer.get_avg_time(), __func__, n, fname, var, seed);
    }
}

void generate_points(vector<Point>& mypoints, size_t totsize, double var, int seed, int blocks, int root, Comm comm)
{
    double t, maxtime, sumtime;
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();
    timer.start_timer();

    BlockPointGenerator<PointTraits> gen(seed, blocks, comm.getcomm());
    gen.generate_points(mypoints, totsize, var);

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] :: []\n", timer.get_max_time(), timer.get_avg_time(), __func__);
    }
}

void write_points_file(const vector<Point>& mypoints, const char *fname, Comm comm)
{
    double t, maxtime, sumtime;
    int myrank = comm.rank();
    int nprocs = comm.size();

    auto timer = comm.get_timer();

    timer.start_timer();
    PointTraits::write_to_file(mypoints.begin(), mypoints.end(), fname, comm);
    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] :: [fname='%s',filesize=%s]\n", timer.get_max_time(), timer.get_avg_time(), __func__, fname, PrettyFileSize::str(fname).c_str());
    }
}
