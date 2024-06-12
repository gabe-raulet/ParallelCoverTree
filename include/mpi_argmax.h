#ifndef MPI_ARGMAX_H_
#define MPI_ARGMAX_H_

#include "mpi_comm.h"

template <class Real, class Index>
struct ArgmaxPair
{
    Index index;
    Real value;

    ArgmaxPair(Index index, Real value) : index(index), value(value) {}

    static void mpi_argmax(void *_in, void *_inout, int *len, MPI_Datatype *dtype)
    {
        ArgmaxPair *in = (ArgmaxPair*)_in;
        ArgmaxPair *inout = (ArgmaxPair*)_inout;

        for (int i = 0; i < *len; ++i)
            if (inout[i].value < in[i].value)
            {
                inout[i].value = in[i].value;
                inout[i].index = in[i].index;
            }
    }

    static void create_mpi_handlers(MPI_Datatype& MPI_ARGMAX_PAIR, MPI_Op& MPI_ARGMAX)
    {
        int blklens[2] = {1,1};
        MPI_Aint disps[2] = {offsetof(ArgmaxPair, index), offsetof(ArgmaxPair, value)};
        MPI_Datatype types[2] = {MPIType<Index>(), MPIType<Real>()};
        MPI_Type_create_struct(2, blklens, disps, types, &MPI_ARGMAX_PAIR);
        MPI_Type_commit(&MPI_ARGMAX_PAIR);
        MPI_Op_create(&mpi_argmax, 1, &MPI_ARGMAX);
    }
};

#endif
