#ifndef TIMERS_H_
#define TIMERS_H_

#include <chrono>
#include <mpi.h>

using namespace std;

struct LocalTimer
{
    using clock = chrono::high_resolution_clock;
    using duration = chrono::duration<double>;
    using time_point = chrono::time_point<clock, duration>;

    time_point start, end;

    static double time_between(const time_point& start, const time_point& end)
    {
        return duration(end - start).count();
    }

    void start_timer() { start = clock::now(); }
    void stop_timer()  { end = clock::now(); }

    double get_elapsed() { return time_between(start, end); }
};

struct MPITimer
{
    int myrank, nprocs;
    MPI_Comm comm;
    double t, maxtime, avgtime;

    MPITimer(MPI_Comm comm) : comm(comm)
    {
        MPI_Comm_rank(comm, &myrank);
        MPI_Comm_size(comm, &nprocs);
    }

    void start_timer()
    {
        MPI_Barrier(comm);
        t = -MPI_Wtime();
    }

    void stop_timer()
    {
        t += MPI_Wtime();
        MPI_Reduce(&t, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&t, &avgtime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avgtime /= nprocs;
    }
};

#endif
