#include <vector>
#include <iostream>
#include <filesystem>
#include <sstream>
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

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed = -MPI_Wtime();

    index_t n = tree.num_points();
    const float *p = tree.getdata();
    int d = tree.getdim();

    index_t chunksize = n / nprocs;
    index_t mychunksize = myrank != nprocs-1? chunksize : n - (nprocs-1)*chunksize;
    index_t myoffset = myrank * chunksize;

    std::vector<std::vector<index_t>> mygraph;
    mygraph.resize(mychunksize);

    const float *myp = &p[d*myoffset];

    for (index_t i = 0; i < mychunksize; ++i)
    {
        mygraph[i] = tree.radii_query(&myp[d*i], epsilon);
    }

    elapsed += MPI_Wtime();

    MPI_Reduce(&elapsed, &maxtime,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed, &proctime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        fprintf(stderr, "Constructed neighborhood graph :: maxtime = %4f (s), avgtime = %4f (s)\n", maxtime, proctime / nprocs);
    }

    if (outfname)
    {
        std::ostringstream ss;
        ss << "rank" << myrank << "_" << outfname;
        std::string fname = ss.str();

        index_t m = 0;

        for (auto it = mygraph.begin(); it != mygraph.end(); ++it)
        {
            m += it->size();
            std::sort(it->begin(), it->end());
        }

        FILE *f = fopen(fname.c_str(), "w");

        fprintf(f, "%lld\t%lld\t%lld\n", n, n, m);

        for (index_t i = 0; i < n; ++i)
        {
            for (index_t j : mygraph[i])
            {
                fprintf(f, "%lld\t%lld\n", i+myoffset, j);
            }
        }

        fclose(f);
    }

    MPI_Finalize();
    return 0;
}
