#include <numeric>

MPICommunicator::MPICommunicator(MPI_Comm comm) : comm(comm)
{
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
}

template <class T>
void MPICommunicator::exscan(const T* sendbuf, T* recvbuf, int count, MPI_Op op)
{
    MPI_Exscan(sendbuf, recvbuf, count, mpi_get_type<T>(), op, comm);
}

template <class T>
void MPICommunicator::exscan(const T& sendbuf, T& recvbuf, MPI_Op op)
{
    exscan(&sendbuf, &recvbuf, 1, op);
}

template <class T>
void MPICommunicator::exscan(const vector<T>& sendbuf, vector<T>& recvbuf, MPI_Op op)
{
    int count = sendbuf.size();
    recvbuf.resize(count);
    exscan(sendbuf.data(), recvbuf.data(), count, op);
}

template <class T>
void MPICommunicator::bcast(T* buffer, int count, int root)
{
    MPI_Bcast(buffer, count, mpi_get_type<T>(), root, comm);
}

template <class T>
void MPICommunicator::bcast(T& buffer, int root)
{
    bcast(&buffer, 1, root);
}

template <class T>
void MPICommunicator::bcast(vector<T>& buffer, int root)
{
    int count = buffer.size();
    bcast(buffer.data(), count, root);
}

template <class T> void MPICommunicator::reduce(const T* sendbuf, T* recvbuf, int count, int root, MPI_Op op)
{
    MPI_Reduce(sendbuf, recvbuf, count, mpi_get_type<T>(), op, root, comm);
}

template <class T> void MPICommunicator::reduce(const T& sendbuf, T& recvbuf, int root, MPI_Op op)
{
    reduce(&sendbuf, &recvbuf, 1, root, op);
}

template <class T> void MPICommunicator::reduce(const vector<T>& sendbuf, vector<T>& recvbuf, int root, MPI_Op op)
{
    int count = sendbuf.size();
    int myrank = getrank();

    if (myrank == root) recvbuf.resize(count);

    reduce(sendbuf.data(), recvbuf.data(), count, root, op);
}

template <class T> void MPICommunicator::allreduce(T* buffer, int count, MPI_Op op)
{
    MPI_Allreduce(MPI_IN_PLACE, buffer, count, mpi_get_type<T>(), op, comm);
}

template <class T> void MPICommunicator::allreduce(T& buffer, MPI_Op op)
{
    allreduce(&buffer, 1, op);
}

template <class T> void MPICommunicator::allreduce(vector<T>& buffer, MPI_Op op)
{
    int count = buffer.size();
    allreduce(buffer.data(), count, op);
}

template <class T> void MPICommunicator::allreduce(const T* sendbuf, T* recvbuf, int count, MPI_Op op)
{
    MPI_Allreduce(sendbuf, recvbuf, count, mpi_get_type<T>(), op, comm);
}

template <class T> void MPICommunicator::allreduce(const T& sendbuf, T* recvbuf, MPI_Op op)
{
    allreduce(&sendbuf, recvbuf, 1, op);
}

template <class T> void MPICommunicator::allreduce(const vector<T>& sendbuf, vector<T>& recvbuf, MPI_Op op)
{
    int count = sendbuf.size();
    recvbuf.resize(count);
    allreduce(sendbuf.data(), recvbuf.data(), count, op);
}

template <class T> void MPICommunicator::gather(const T* sendbuf, int count, T* recvbuf, int root)
{
    MPI_Gather(sendbuf, count, mpi_get_type<T>(), recvbuf, count, mpi_get_type<T>(), root, comm);
}

template <class T> void MPICommunicator::gather(const T& sendbuf, T* recvbuf, int root)
{
    gather(&sendbuf, 1, recvbuf, root);
}

template <class T> void MPICommunicator::gather(const vector<T>& sendbuf, vector<T>& recvbuf, int root)
{
    int count = sendbuf.size();
    if (myrank == root) recvbuf.resize(count * nprocs);
    gather(sendbuf.data(), count, recvbuf.data(), root);
}

template <class T> void MPICommunicator::gatherv(const vector<T>& sendbuf, vector<T>& recvbuf, int root)
{
    vector<int> recvcounts, displs;
    int sendcount = sendbuf.size();

    gather({sendcount}, recvcounts, root);

    if (myrank == root)
    {
        displs.resize(nprocs);
        displs.front() = 0;
        partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
        recvbuf.resize(recvcounts.back() + displs.back());
    }

    MPI_Gatherv(sendbuf.data(), sendcount, mpi_get_type<T>(), recvbuf.data(), recvcounts.data(), displs.data(), mpi_get_type<T>(), root, comm);
}

template <class T> void MPICommunicator::allgather(T* buffer, int count)
{
    MPI_Allgather(MPI_IN_PLACE, count, mpi_get_type<T>(), buffer, count, mpi_get_type<T>(), comm);
}

template <class T> void MPICommunicator::allgather(T* buffer)
{
    allgather(buffer, 1);
}

template <class T> void MPICommunicator::allgather(vector<T>& buffer)
{
    int count = buffer.size();
    allgather(buffer.data(), count);
}

template <class T> void MPICommunicator::allgather(const T* sendbuf, int count, T* recvbuf)
{
    MPI_Allgather(sendbuf, count, mpi_get_type<T>(), recvbuf, count, mpi_get_type<T>(), comm);
}

template <class T> void MPICommunicator::allgather(const T& sendbuf, T* recvbuf)
{
    allgather(&sendbuf, 1, recvbuf);
}

template <class T> void MPICommunicator::allgather(const T& sendbuf, vector<T>& recvbuf)
{
    recvbuf.resize(nprocs);
    allgather(sendbuf, recvbuf.data());
}

template <class T> void MPICommunicator::allgatherv(const vector<T>& sendbuf, vector<T>& recvbuf)
{
    int sendcount = sendbuf.size();
    vector<int> recvcounts(nprocs), displs(nprocs);

    allgather(sendcount, recvcounts);

    displs.front() = 0;
    partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
    recvbuf.resize(recvcounts.back() + displs.back());

    MPI_Allgatherv(sendbuf.data(), sendcount, mpi_get_type<T>(), recvbuf.data(), recvcounts.data(), displs.data(), mpi_get_type<T>(), comm);
}

template <class T> void MPICommunicator::scatter(const T* sendbuf, int count, T* recvbuf, int root)
{
    MPI_Scatter(sendbuf, count, mpi_get_type<T>(), recvbuf, count, mpi_get_type<T>(), root, comm);
}

template <class T> void MPICommunicator::scatter(const T* sendbuf, T& recvbuf, int root)
{
    scatter(sendbuf, 1, &recvbuf, root);
}

template <class T> void MPICommunicator::scatter(const vector<T>& sendbuf, vector<T>& recvbuf, int root)
{
    int count = sendbuf.size();
    if (myrank == root) recvbuf.resize(count * nprocs);
    scatter(sendbuf.data(), count, recvbuf.data(), root);
}

template <class T> void MPICommunicator::scatterv(const vector<T>& sendbuf, const vector<int>& sendcounts, vector<T>& recvbuf, int root)
{
    int recvcount;
    scatter(sendcounts.data(), recvcount, root);
    recvbuf.resize(recvcount);

    vector<int> displs;
    if (myrank == root)
    {
        displs.resize(nprocs);
        displs.front() = 0;
        partial_sum(sendcounts.begin(), sendcounts.end()-1, displs.begin()+1);
    }

    MPI_Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), mpi_get_type<T>(), recvbuf.data(), recvcount, mpi_get_type<T>(), root, comm);
}

bool MPICommunicator::is_self_comm() const
{
    int result;
    MPI_Comm_compare(comm, MPI_COMM_SELF, &result);
    return static_cast<bool>(result == MPI_IDENT);
}
