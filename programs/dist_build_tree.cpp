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
#include <mpi.h>
#include "DistCoverTree.h"
#include "VectorIO.h"
#include "read_args.h"

int main(int argc, char *argv[])
{
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    double base = 2.;
    char *infname = NULL, *outfname = NULL;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        if (myrank == 0)
        {
            fprintf(stderr, "Usage: %s [options] [..]\n", argv[0]);
            fprintf(stderr, "Options: -i FILE    input vector filename [required]\n");
            fprintf(stderr, "         -o FILE    output cover tree filename [required]\n");
            fprintf(stderr, "         -b FLOAT   cover tree base [default: %.2f]\n", base);
            fprintf(stderr, "         -h         help message\n");
        }

        MPI_Finalize();
        return 0;
    }

    infname = read_string_arg(argc, argv, "-i", NULL);
    outfname = read_string_arg(argc, argv, "-o", NULL);
    base = read_double_arg(argc, argv, "-b", &base);

    double elapsed, maxtime, proctime;

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed = -MPI_Wtime();

    int d;
    index_t n;
    std::vector<float> pointmem;

    if (myrank == 0)
    {
        size_t _n;
        pointmem = read_vecs_file(infname, &d, &_n);
        n = _n;
    }

    MPI_Bcast(&d, 1, MPI_INT,     0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);

    if (myrank != 0)
    {
        pointmem.resize(n*d);
    }

    MPI_Bcast(pointmem.data(), static_cast<int>(n*d), MPI_FLOAT, 0, MPI_COMM_WORLD);

    elapsed += MPI_Wtime();

    MPI_Reduce(&elapsed, &maxtime,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myrank == 0) fprintf(stderr, "Read points file '%s' (%lld points of dimension %d) :: maxtime = %.4f (s), avgtime = %.4f (s)\n", infname, n, d, maxtime, proctime / nprocs);

    MPI_Finalize();
    return 0;
}
