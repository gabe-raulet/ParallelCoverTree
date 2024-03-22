#include "MPITimer.h"

MPITimer::MPITimer(MPI_Comm comm, int root) : comm(comm), root(root)
{
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
}

void MPITimer::start_timer()
{
    MPI_Barrier(comm);
    t = -MPI_Wtime();
}

void MPITimer::stop_timer()
{
    t += MPI_Wtime();
    MPI_Reduce(&t, &maxtime,  1, MPI_DOUBLE, MPI_MAX, root, comm);
    MPI_Reduce(&t, &proctime, 1, MPI_DOUBLE, MPI_SUM, root, comm);
}
