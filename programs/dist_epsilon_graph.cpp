#include "DistCoverTree.h"
#include "Point.h"
#include "MPITimer.h"
#include "read_args.h"
#include "version.h"
#include <iostream>
#include <sstream>
#include <chrono>
#include <assert.h>
#include <stdio.h>
#include <mpi.h>

using namespace std;

string return_current_time_and_date();
string program_str(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int64_t n;
    double var = 10.0;
    double base = 2.0;
    int seed = -1;
    int verbose = 0;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        if (!myrank)
        {
            fprintf(stderr, "Usage: %s [options]\n", argv[0]);
            fprintf(stderr, "Options: -n INT    number of points [required]\n");
            fprintf(stderr, "         -V FLOAT  variance [default: %.2f]\n", var);
            fprintf(stderr, "         -C FLOAT  cover base [default: %.2f]\n", base);
            fprintf(stderr, "         -s INT    rng seed [default: random]\n");
            fprintf(stderr, "         -v        verbose\n");
            fprintf(stderr, "         -h        help message\n");
        }

        MPI_Finalize();
        return 0;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    var = read_double_arg(argc, argv, "-V", &var);
    seed = read_int_arg(argc, argv, "-s", &seed);
    base = read_double_arg(argc, argv, "-C", &base);
    verbose = !!(find_arg_idx(argc, argv, "-v") >= 0);

    if (seed < 0)
    {
        random_device rd;
        default_random_engine gen(rd());
        uniform_int_distribution<int> dis{0, numeric_limits<int>::max()};
        seed = dis(gen);
    }

    if (!myrank) fprintf(stderr, "program: %s\ncommit: " GIT_COMMIT "\ndate and time: %s\nMPI ranks: %d\n\n", program_str(argc, argv).c_str(), return_current_time_and_date().c_str(), nprocs);

    MPITimer timer(MPI_COMM_WORLD, 0);

    timer.start_timer();

    vector<Point> mypoints = Point::dist_random_points(n, var, seed, 0, MPI_COMM_WORLD);

    timer.stop_timer();

    if (!myrank) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (dist_random_points) :: [n=%lld,var=%.2f,seed=%d]\n", timer.get_max_time(), timer.get_avg_time(), n, var, seed);

    timer.start_timer();

    DistCoverTree tree(mypoints, base, MPI_COMM_WORLD);
    tree.build_tree(static_cast<bool>(verbose));

    timer.stop_timer();
    tree.print_timing_results();

    if (!myrank) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (dist_build_tree) [num_vertices=%lld,num_levels=%lld,base=%.2f]\n", timer.get_max_time(), timer.get_avg_time(), tree.num_vertices(), tree.num_levels(), base);

    MPI_Finalize();
    return 0;
}

string program_str(int argc, char *argv[])
{
    stringstream ss;

    for (int i = 0; i < argc; ++i)
        ss << argv[i] << " ";

    return ss.str();
}

string return_current_time_and_date()
{
    /*
     * Shamelessly copied from here: https://stackoverflow.com/questions/17223096/outputting-date-and-time-in-c-using-stdchrono
     */

    auto now = chrono::system_clock::now();
    auto in_time_t = chrono::system_clock::to_time_t(now);

    stringstream ss;
    ss << put_time(localtime(&in_time_t), "%Y/%m/%d %X");
    return ss.str();
}
