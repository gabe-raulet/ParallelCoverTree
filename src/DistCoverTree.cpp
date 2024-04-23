#include "DistCoverTree.h"
#include "Point.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <limits>
#include <iomanip>
#include <cassert>
#include <stdio.h>

DistCoverTree::DistCoverTree(const vector<Point>& mypoints, double base, MPI_Comm comm)
    : max_radius(-1), base(base), mysize(mypoints.size()), mypoints(mypoints), comm(comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Exscan(&mysize, &myoffset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (myrank == 0) myoffset = 0;

    totsize = mysize + myoffset;
    MPI_Bcast(&totsize, 1, MPI_INT64_T, nprocs-1, comm);
}

void DistCoverTree::build_tree() {}

int64_t DistCoverTree::add_vertex(int64_t point_id, int64_t parent_id)
{
    int64_t vertex_level;
    int64_t vertex_id = pt.size(); // this is probably wrong

    pt.push_back(point_id);
    children.emplace_back();

    if (parent_id >= 0)
    {
        vertex_level = level[parent_id] + 1;
        children[parent_id].push_back(vertex_id);
    }
    else vertex_level = 0;

    level.push_back(vertex_level);
    return vertex_id;
}

double DistCoverTree::vertex_ball_radius(int64_t vertex_id) const
{
    return pow(base, -1. * level[vertex_id]);
}

void DistCoverTree::initialize_root_hub()
{
    my_dists.resize(mysize);
    my_hub_vtx_ids.resize(mysize);
    my_hub_pt_ids.resize(mysize);

    Point root_pt = mypoints.front();

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);
    MPI_Type_commit(&MPI_POINT);
    MPI_Bcast(&root_pt, 1, MPI_POINT, 0, comm);
    MPI_Type_free(&MPI_POINT);

    for (int64_t i = 0; i < mysize; ++i)
    {
        my_dists[i] = root_pt.distance(mypoints[i]);
        my_hub_vtx_ids[i] = my_hub_pt_ids[i] = 0;
        max_radius = max(my_dists[i], max_radius);
    }

    MPI_Allreduce(MPI_IN_PLACE, &max_radius, 1, MPI_DOUBLE, MPI_MAX, comm);
}

struct ArgmaxPair
{
    int64_t index;
    double value;

    ArgmaxPair(int64_t index, double value) : index(index), value(value) {}

    static void create_mpi_dtype(MPI_Datatype *dtype)
    {
        int blklens[2] = {1,1};
        MPI_Aint disps[2] = {offsetof(ArgmaxPair, index), offsetof(ArgmaxPair, value)};
        MPI_Datatype types[2] = {MPI_INT64_T, MPI_DOUBLE};
        MPI_Type_create_struct(2, blklens, disps, types, dtype);
        MPI_Type_commit(dtype);
    }

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

    static void create_mpi_op(MPI_Op *op)
    {
        MPI_Op_create(&mpi_argmax, 1, op);
    }
};


void DistCoverTree::compute_farthest_hub_pts()
{
    unordered_map<int64_t, pair<int64_t, double>> my_argmaxes; // maps hub id to (point id, distance) pair
    transform(hub_chains.begin(), hub_chains.end(), inserter(my_argmaxes, my_argmaxes.end()),
            [](const auto& chain) { return make_pair(chain.first, make_pair(-1, -1.0)); });

    for (int64_t i = 0; i < mysize; ++i)
    {
        int64_t hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            auto& it = my_argmaxes.find(hub_id)->second;

            if (my_dists[i] > it.second)
            {
                it.first = i;
                it.second = my_dists[i];
            }
        }
    }

    vector<int64_t> hub_ids;
    vector<ArgmaxPair> my_argmax_pairs;

    hub_ids.reserve(my_argmaxes.size());
    my_argmax_pairs.reserve(my_argmaxes.size());

    for (auto it = my_argmaxes.begin(); it != my_argmaxes.end(); ++it)
    {
        hub_ids.push_back(it->first);
        my_argmax_pairs.emplace_back(it->second.first + myoffset, it->second.second);
    }

    MPI_Op MPI_ARGMAX;
    MPI_Datatype MPI_ARGMAX_PAIR;

    ArgmaxPair::create_mpi_op(&MPI_ARGMAX);
    ArgmaxPair::create_mpi_dtype(&MPI_ARGMAX_PAIR);

    MPI_Allreduce(MPI_IN_PLACE, my_argmax_pairs.data(), static_cast<int>(my_argmax_pairs.size()), MPI_ARGMAX_PAIR, MPI_ARGMAX, comm);

    MPI_Op_free(&MPI_ARGMAX);
    MPI_Type_free(&MPI_ARGMAX_PAIR);

    farthest_hub_pts.clear();

    for (int64_t i = 0; i < my_argmax_pairs.size(); ++i)
    {
        farthest_hub_pts.insert({hub_ids[i], {my_argmax_pairs[i].index, my_argmax_pairs[i].value}});
    }
}
