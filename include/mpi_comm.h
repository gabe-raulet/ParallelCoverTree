#ifndef MPI_COMM_H_
#define MPI_COMM_H_

#include <vector>
#include <map>
#include <cstdint>
#include <cassert>
#include <mpi.h>

using namespace std;

struct type_info_compare
{
    bool operator()(type_info const* lhs, type_info const* rhs) const
    {
        return lhs->before(*rhs);
    }
};

class MPIDatatypeCache
{
    private:

        typedef map<type_info const*, MPI_Datatype, type_info_compare> stored_map_type;
        stored_map_type map;

    public:

        void clear()
        {
            int is_finalized = 0;
            MPI_Finalized(&is_finalized);

            if (!is_finalized)
                for (auto it = map.begin(); it != map.end(); ++it)
                    MPI_Type_free(&(it->second));
        }

        ~MPIDatatypeCache() { clear(); }

        MPI_Datatype get(const type_info* t)
        {
            auto pos = map.find(t);
            return pos != map.end()? pos->second : MPI_DATATYPE_NULL;
        }

        void set(const type_info* t, MPI_Datatype datatype)
        {
            map[t] = datatype;
        }
};

MPIDatatypeCache mpidtc;

template <class T>
MPI_Datatype MPIType()
{
    type_info const* t = &typeid(T);
    MPI_Datatype datatype = mpidtc.get(t);

    if (datatype == MPI_DATATYPE_NULL)
    {
        MPI_Type_contiguous(sizeof(T), MPI_CHAR, &datatype);
        MPI_Type_commit(&datatype);
        mpidtc.set(t, datatype);
    }

    return datatype;
}

template <> MPI_Datatype MPIType<char>() { return MPI_CHAR; }
template <> MPI_Datatype MPIType<signed char>() { return MPI_SIGNED_CHAR; }
template <> MPI_Datatype MPIType<short>() { return MPI_SHORT; }
template <> MPI_Datatype MPIType<int>() { return MPI_INT; }
template <> MPI_Datatype MPIType<long>() { return MPI_LONG; }
template <> MPI_Datatype MPIType<long long>() { return MPI_LONG_LONG; }
template <> MPI_Datatype MPIType<unsigned char>() { return MPI_UNSIGNED_CHAR; }
template <> MPI_Datatype MPIType<unsigned short>() { return MPI_UNSIGNED_SHORT; }
template <> MPI_Datatype MPIType<unsigned int>() { return MPI_UNSIGNED; }
template <> MPI_Datatype MPIType<unsigned long>() { return MPI_UNSIGNED_LONG; }
template <> MPI_Datatype MPIType<unsigned long long>() { return MPI_UNSIGNED_LONG_LONG; }
template <> MPI_Datatype MPIType<float>() { return MPI_FLOAT; }
template <> MPI_Datatype MPIType<double>() { return MPI_DOUBLE; }
template <> MPI_Datatype MPIType<long double>() { return MPI_LONG_DOUBLE; }

class MPICommunicator
{
    public:

        MPICommunicator(MPI_Comm comm);

        bool is_self_comm() const;

        MPI_Comm getcomm() const { return comm; }
        int getrank() const { return myrank; }
        int getsize() const { return nprocs; }

        void barrier() { MPI_Barrier(comm); }

        template <class T> void exscan(const T& sendbuf, T& recvbuf, MPI_Op op);
        template <class T> void exscan(const T* sendbuf, T* recvbuf, int count, MPI_Op op);
        template <class T> void exscan(const vector<T>& sendbuf, vector<T>& recvbuf, MPI_Op op);

        template <class T> void bcast(T& buffer, int root);
        template <class T> void bcast(T* buffer, int count, int root);
        template <class T> void bcast(vector<T>& buffer, int root);

        template <class T> void reduce(const T* sendbuf, T* recvbuf, int count, int root, MPI_Op op);
        template <class T> void reduce(const T& sendbuf, T& recvbuf, int root, MPI_Op op);
        template <class T> void reduce(const vector<T>& sendbuf, vector<T>& recvbuf, int root, MPI_Op op);

        template <class T> void allreduce(T& buffer, MPI_Op op);
        template <class T> void allreduce(T* buffer, int count, MPI_Op op);
        template <class T> void allreduce(vector<T>& buffer, MPI_Op op);
        template <class T> void allreduce(const T& sendbuf, T* recvbuf, MPI_Op op);
        template <class T> void allreduce(const T* sendbuf, T* recvbuf, int count, MPI_Op op);
        template <class T> void allreduce(const vector<T>& sendbuf, vector<T>& recvbuf, MPI_Op op);

        template <class T> void gather(const T& sendbuf, T* recvbuf, int root);
        template <class T> void gather(const T* sendbuf, int sendcount, T* recvbuf, int root);
        template <class T> void gather(const vector<T>& sendbuf, vector<T>& recvbuf, int root);
        template <class T> void gatherv(const vector<T>& sendbuf, vector<T>& recvbuf, int root);

        template <class T> void allgather(T* buffer);
        template <class T> void allgather(T* buffer, int sendcount);
        template <class T> void allgather(vector<T>& buffer);
        template <class T> void allgather(const T& sendbuf, T* recvbuf);
        template <class T> void allgather(const T& sendbuf, vector<T>& recvbuf);
        template <class T> void allgather(const T* sendbuf, int sendcount, T* recvbuf);
        template <class T> void allgatherv(const vector<T>& sendbuf, vector<T>& recvbuf);

        template <class T> void scatter(const T* sendbuf, T& recvbuf, int root);
        template <class T> void scatter(const T* sendbuf, int sendcount, T* recvbuf, int root);
        template <class T> void scatter(const vector<T>& sendbuf, vector<T>& recvbuf, int root);
        template <class T> void scatterv(const vector<T>& sendbuf, const vector<int>& sendcounts, vector<T>& recvbuf, int root);

        template <class T>
        static MPI_Datatype mpi_get_type() { return MPIType<T>(); }

    private:

        int myrank, nprocs;
        MPI_Comm comm;
};

#include "mpi_comm.hpp"

#endif
