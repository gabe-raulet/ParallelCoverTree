#ifndef MPI_TIMER_H_
#define MPI_TIMER_H_

#include <mpi.h>

struct MPITimer
{
    MPITimer(MPI_Comm comm, int root);

    void start_timer();
    void stop_timer();

    double get_max_time() const { return maxtime; }
    double get_proc_time() const { return proctime; }
    double get_avg_time() const { return proctime / nprocs; }

    MPI_Comm comm;
    int root, myrank, nprocs;
    double t, maxtime, proctime;
};

#endif
