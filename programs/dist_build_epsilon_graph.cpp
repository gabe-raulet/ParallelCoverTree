#include <vector>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <sys/stat.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <mpi.h>
#include "CoverTree.h"
#include "VectorIO.h"
#include "read_args.h"

int main(int argc, char *argv[])
{
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    double epsilon;
    char *infname = NULL;
    char *outfname = NULL;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        if (myrank == 0)
        {
            fprintf(stderr, "Usage: %s [options] [..]\n", argv[0]);
            fprintf(stderr, "Options: -i FILE    input cover tree filename [required]\n");
            fprintf(stderr, "         -r FLOAT   epsilon radius [required]\n");
            fprintf(stderr, "         -o FILE    neighborhood graph file [default: none]\n");
            fprintf(stderr, "         -h         help message\n");
        }

        MPI_Finalize();
        return 0;
    }

    infname = read_string_arg(argc, argv, "-i", NULL);
    epsilon = read_double_arg(argc, argv, "-r", NULL);
    outfname = read_string_arg(argc, argv, "-o", &outfname);

    double elapsed, maxtime, proctime;

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed = -MPI_Wtime();

    CoverTree tree;
    tree.read_from_file(infname);

    elapsed += MPI_Wtime();

    MPI_Reduce(&elapsed, &maxtime,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        fprintf(stderr, "Read cover tree from file '%s' :: maxtime = %4f (s), avgtime = %4f (s)\n", infname, maxtime, proctime / nprocs);
    }

    MPI_Finalize();
    return 0;
}
