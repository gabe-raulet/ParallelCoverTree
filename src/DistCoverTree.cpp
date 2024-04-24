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

void DistCoverTree::build_tree()
{
    initialize_root_hub();

    while (!hub_chains.empty())
    {
        compute_farthest_hub_pts();
        update_hub_chains();
        process_leaf_chains();
        process_split_chains();
        update_dists_and_pointers();
    }
}

int64_t DistCoverTree::add_vertex(int64_t point_id, int64_t parent_id)
{
    int64_t vertex_level;
    int64_t vertex_id = pt.size();

    pt.push_back(point_id);
    children.emplace_back();

    if (parent_id >= 0)
    {
        vertex_level = level[parent_id] + 1;
        children[parent_id].push_back(vertex_id);
    }
    else vertex_level = 0;

    level.push_back(vertex_level);
    cover_map.insert({vertex_id, pow(base, -1. * vertex_level)});

    return vertex_id;
}

double DistCoverTree::vertex_ball_radius(int64_t vertex_id) const
{
    return cover_map.find(vertex_id)->second;
}

void DistCoverTree::initialize_root_hub()
{
    my_dists.resize(mysize);
    my_hub_vtx_ids.resize(mysize);
    my_hub_pt_ids.resize(mysize);

    Point root_pt = mypoints.front();
    int64_t root_id = add_vertex(0, -1);

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);
    MPI_Bcast(&root_pt, 1, MPI_POINT, 0, comm);
    MPI_Type_free(&MPI_POINT);

    for (int64_t i = 0; i < mysize; ++i)
    {
        my_dists[i] = root_pt.distance(mypoints[i]);
        my_hub_vtx_ids[i] = root_id;
        my_hub_pt_ids[i] = pt[root_id];
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

void DistCoverTree::update_hub_chains()
{
    int64_t hub_id;
    pair<int64_t, double> farthest_pt;
    split_chains.clear(), leaf_chains.clear();

    for (auto it = farthest_hub_pts.begin(); it != farthest_hub_pts.end(); ++it)
    {
        tie(hub_id, farthest_pt) = *it;
        int64_t farthest_pt_id = farthest_pt.first;
        double farthest_dist = farthest_pt.second / max_radius;

        if (farthest_dist == 0)
        {
            hub_chains.erase(hub_id);
            leaf_chains.insert(hub_id);
        }
        else if (farthest_dist <= (vertex_ball_radius(hub_id) / base))
        {
            split_chains.push_back(hub_id);
        }
        else
        {
            hub_chains.find(hub_id)->second.push_back(farthest_pt_id);
        }
    }
}

int64_t DistCoverTree::batch_new_vertex(int64_t point_id, int64_t parent_id)
{
    my_new_vertex_pt_ids.push_back(point_id);
    my_new_vertex_hub_ids.push_back(parent_id);
    return pt.size() + my_new_vertex_pt_ids.size();
}

void DistCoverTree::add_batched_vertices()
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    vector<int64_t> new_vertex_pt_ids, new_vertex_hub_ids;

    vector<int> recvcounts(nprocs);
    vector<int> displs(nprocs);
    int sendcount = static_cast<int>(my_new_vertex_pt_ids.size());

    recvcounts[myrank] = sendcount;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    displs.front() = 0;
    partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    int totrecv = recvcounts.back() + displs.back();
    new_vertex_pt_ids.resize(totrecv);
    new_vertex_hub_ids.resize(totrecv);

    MPI_Allgatherv(my_new_vertex_pt_ids.data(), sendcount, MPI_INT64_T, new_vertex_pt_ids.data(), recvcounts.data(), displs.data(), MPI_INT64_T, comm);
    MPI_Allgatherv(my_new_vertex_hub_ids.data(), sendcount, MPI_INT64_T, new_vertex_hub_ids.data(), recvcounts.data(), displs.data(), MPI_INT64_T, comm);

    for (int i = 0; i < totrecv; ++i)
    {
        add_vertex(new_vertex_pt_ids[i], new_vertex_hub_ids[i]);
    }

    my_new_vertex_pt_ids.clear();
    my_new_vertex_hub_ids.clear();
}

void DistCoverTree::process_leaf_chains()
{
    if (!leaf_chains.empty())
    {
        for (int64_t i = 0; i < mysize; ++i)
        {
            int64_t hub_id = my_hub_vtx_ids[i];

            if (leaf_chains.find(hub_id) != leaf_chains.end())
            {
                batch_new_vertex(i + myoffset, hub_id);
                my_hub_vtx_ids[i] = my_hub_pt_ids[i] = -1;
                my_dists[i] = 0;
            }
        }
    }

    add_batched_vertices();
}

void DistCoverTree::process_split_chains()
{
    if (!split_chains.empty())
    {
        unordered_map<int64_t, int64_t> hub_pt_id_updates;
        vector<int64_t> new_hub_pts, new_hub_ids;
        vector<vector<int64_t>*> chain_ptrs(split_chains.size(), nullptr);
        transform(split_chains.begin(), split_chains.end(), chain_ptrs.begin(),
                [&](int64_t hub_id) { return &hub_chains.find(hub_id)->second; });

        int64_t num_new_hub_pts = 0;
        for (auto& chain_ptr : chain_ptrs)
            num_new_hub_pts += chain_ptr->size();

        new_hub_pts.reserve(num_new_hub_pts);
        new_hub_ids.reserve(num_new_hub_pts);

        for (int64_t i = 0; i < chain_ptrs.size(); ++i)
        {
            new_hub_pts.insert(new_hub_pts.end(), chain_ptrs[i]->cbegin(), chain_ptrs[i]->cend());
            new_hub_ids.insert(new_hub_ids.end(), chain_ptrs[i]->size(), split_chains[i]);
            hub_chains.erase(split_chains[i]);
        }

        for (int64_t i = 0; i < new_hub_pts.size(); ++i)
        {
            int64_t vtx_id = add_vertex(new_hub_pts[i], new_hub_ids[i]);
            hub_chains.insert({vtx_id, {new_hub_pts[i]}});
            hub_pt_id_updates.insert({new_hub_pts[i], vtx_id});
        }

        for (int64_t i = 0; i < mysize; ++i)
        {
            int64_t closest_pt_id = my_hub_pt_ids[i];
            auto it = hub_pt_id_updates.find(closest_pt_id);

            if (it != hub_pt_id_updates.end())
                my_hub_vtx_ids[i] = it->second;
        }
    }
}

unordered_map<int64_t, Point> DistCoverTree::collect_points(const vector<int64_t>& point_ids) const
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int sendcount;
    vector<int> recvcounts(nprocs), displs(nprocs);
    vector<Point> sendbuf, recvbuf;
    vector<int64_t> myids, ids;

    for (int64_t id : point_ids)
        if (myoffset <= id && id < myoffset + mysize)
        {
            sendbuf.push_back(mypoints[id]);
            myids.push_back(id);
        }

    sendcount = static_cast<int>(sendbuf.size());

    MPI_Allgather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    displs.front() = 0;
    partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    int totrecv = recvcounts.back() + displs.back();
    recvbuf.resize(totrecv);
    ids.resize(totrecv);

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);

    MPI_Allgatherv(sendbuf.data(), sendcount, MPI_POINT, recvbuf.data(), recvcounts.data(), displs.data(), MPI_POINT, comm);
    MPI_Allgatherv(myids.data(), sendcount, MPI_INT64_T, ids.data(), recvcounts.data(), displs.data(), MPI_INT64_T, comm);

    unordered_map<int64_t, Point> collected_points;
    collected_points.reserve(totrecv);

    for (int64_t i = 0; i < static_cast<int64_t>(totrecv); ++i)
        collected_points.insert({ids[i], recvbuf[i]});

    return collected_points;
}

void DistCoverTree::update_dists_and_pointers()
{
    vector<int64_t> point_ids;

    for (int64_t i = 0; i < mysize; ++i)
    {
        int64_t hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            int64_t last_chain_pt_id = hub_chains.find(hub_id)->second.back();
            point_ids.push_back(last_chain_pt_id);
        }
    }

    unordered_map<int64_t, Point> collected_points = collect_points(point_ids);

    for (int64_t i = 0; i < mysize; ++i)
    {
        int64_t hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            int64_t last_chain_pt_id = hub_chains.find(hub_id)->second.back();
            double lastdist = my_dists[i];
            double curdist = mypoints[i].distance(collected_points.find(last_chain_pt_id)->second);

            if (curdist <= lastdist)
            {
                my_dists[i] = curdist;
                my_hub_pt_ids[i] = last_chain_pt_id;
            }
        }
    }
}
