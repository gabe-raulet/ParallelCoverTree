#ifndef MPI_ENV_H_
#define MPI_ENV_H_

#include <map>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <string>
#include <sstream>
#include <numeric>
#include <mpi.h>
#include "mytypeinfo.h"
#include "misc.h"

namespace MPIEnv
{
    using namespace std;

    int initialize(int *argc, char **argv[]);
    int finalize();
    void exit(int err);

    bool is_initialized();
    bool is_finalized();

    bool comms_equal(MPI_Comm lhs, MPI_Comm rhs);

    struct type_info_compare
    {
        bool operator()(const type_info *lhs, const type_info *rhs) const;
    };

    class TypeCache
    {
        public:

            TypeCache();
            ~TypeCache();

            MPI_Datatype get_type(const type_info *t);
            void set_type(const type_info *t, MPI_Datatype dtype);
            void free_types();
            void reset();

        private:

            map<const type_info*, MPI_Datatype, type_info_compare> type_map;
    };

    class CommTimer
    {
        public:

            CommTimer(MPI_Comm world)
            {
                MPI_Comm_dup(world, &comm);
                MPI_Comm_rank(comm, &myrank);
                MPI_Comm_size(comm, &nprocs);
            }

            CommTimer(const CommTimer& rhs)
            {
                MPI_Comm_dup(rhs.comm, &comm);
                MPI_Comm_rank(comm, &myrank);
                MPI_Comm_size(comm, &nprocs);
            }

            ~CommTimer() { MPI_Comm_free(&comm); }

            void start_timer()
            {
                MPI_Barrier(comm);
                t = -MPI_Wtime();
            }

            void stop_timer()
            {
                t += MPI_Wtime();
                MPI_Reduce(&t, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
                MPI_Reduce(&t, &sumtime, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
            }

            double get_my_time() const { return t; }
            double get_max_time() const { return maxtime; }
            double get_sum_time() const { return sumtime; }
            double get_avg_time() const { return sumtime/nprocs; }

        private:

            int myrank, nprocs;
            MPI_Comm comm;
            double t, maxtime, sumtime;
    };

    class Comm
    {
        public:

            static Comm comm_world() { return Comm(MPI_COMM_WORLD); }
            static Comm comm_self() { return Comm(MPI_COMM_SELF); }
            static Comm comm_null() { return Comm(MPI_COMM_NULL); }

            bool operator==(const Comm& rhs) const;
            bool operator!=(const Comm& rhs) const { return !(*this == rhs); }

            Comm& operator=(const Comm& rhs);

            MPI_Comm getcomm() const { return commbuf[0]; }
            int rank() const { return myrank; }
            int size() const { return nprocs; }
            tuple <int, int, MPI_Comm> comminfo() const { return {myrank, nprocs, commbuf[0]}; }

            CommTimer get_timer() const { return CommTimer(getcomm()); }

            bool is_distributed() const { return (nprocs > 1); }

            Comm();
            Comm(MPI_Comm comm);
            Comm(const Comm& rhs);
            ~Comm();

            void log_at_root_rank(ostream& os, const char *msg, const char *func) const;
            void log_at_root_rank(ostream& os, const char *msg, const char *func, const CommTimer& timer) const;

            void log_at_all_ranks(ostream& os, const char *msg, const char *func) const;
            void log_at_all_ranks(ostream& os, const char *msg, const char *func, const CommTimer& timer) const;

            void swap(Comm& rhs) noexcept;

            template <class T>
            bool is_same_val(T val) const;

            template <class T>
            bool are_same_vals(const vector<T>& vals) const;

            int barrier() const;

            template <class T> int reduce(const T* sendbuf, T* recvbuf, int count, int root, MPI_Op op) const;
            template <class T> int reduce(const T& sendbuf, T& recvbuf, int root, MPI_Op op) const;
            template <class T> int reduce(const vector<T>& sendbuf, vector<T>& recvbuf, int root, MPI_Op op) const;
            template <class T> int reduce(const vector<T>& sendbuf, vector<T>& recvbuf, int root, const vector<MPI_Op>& ops) const;

            template <class T> int bcast(T* buffer, int count, int root) const;
            template <class T> int bcast(T& buffer, int root) const;
            template <class T> int bcast(vector<T>& buffer, int root) const;

            template <class T> int exscan(const T* sendbuf, T* recvbuf, int count, MPI_Op op, T identity) const;
            template <class T> int exscan(const T& sendbuf, T& recvbuf, MPI_Op op, T identity) const;

            template <class T> int allreduce(const T* sendbuf, T* recvbuf, int count, MPI_Op op) const;
            template <class T> int allreduce(const T& sendbuf, T& recvbuf, MPI_Op op) const;

            template <class T> int allreduce(T* buffer, int count, MPI_Op op) const;
            template <class T> int allreduce(T& buffer, MPI_Op op) const;

            template <class T> int gatherv(const vector<T>& sendbuf, vector<T>& recvbuf, int root) const;
            template <class T> int allgatherv(const vector<T>& sendbuf, vector<T>& recvbuf) const;

            template <class T> int scatterv(const vector<T>& sendbuf, const vector<int>& sendcounts, vector<T>& recvbuf, int root) const;

            template <string_type String>
            int file_write_at_all_once(const char *fname, const String& mybuf, bool truncate=true) const;
            int file_write_at_all_once(const char *fname, const char *mybuf, int count, bool truncate=true) const;

            template <string_type String>
            int file_read_at_all_once(const char *fname, String& mybuf, int unitsize) const;

        private:

            int myrank, nprocs;
            MPI_Comm commbuf[1];

            void init(MPI_Comm comm);
    };

    template <class Real, class Index>
    struct ArgmaxPair
    {
        Index index;
        Real value;

        ArgmaxPair(Index index, Real value) : index(index), value(value) {}

        static void mpi_argmax(void *_in, void *_inout, int *len, MPI_Datatype *dtype);
        static void create_mpi_handlers(MPI_Datatype& MPI_ARGMAX_PAIR, MPI_Op& MPI_ARGMAX);
    };

    TypeCache cache;

    template <class T> MPI_Datatype mpi_type();
    template <class T> MPI_Datatype mpi_commit_type();
};

#include "mpi_env.hpp"

#endif
