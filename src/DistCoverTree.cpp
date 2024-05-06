#include "DistCoverTree.h"
#include "Point.h"
#include "MPITimer.h"
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
    : max_radius(-1),
      base(base),
      mysize(mypoints.size()),
      mypoints(mypoints),
      niters(0),
      nlevels(0),
      comm(comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Exscan(&mysize, &myoffset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (myrank == 0) myoffset = 0;

    totsize = mysize + myoffset;
    MPI_Bcast(&totsize, 1, MPI_INT64_T, nprocs-1, comm);
}

void DistCoverTree::set_times_to_zero()
{
    overall_time = 0;
    initialize_root_hub_time = 0;
    compute_farthest_hub_pts_time = 0;
    update_hub_chains_time = 0;
    process_leaf_chains_time = 0;
    process_split_chains_time = 0;
    update_dists_and_pointers_time = 0;
}

unordered_map<int64_t, int64_t> DistCoverTree::get_hub_counts() const
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    unordered_map<int64_t, size_t> hub_idmap;

    for (auto it = hub_chains.begin(); it != hub_chains.end(); ++it)
    {
        hub_idmap.insert({it->first, hub_idmap.size()});
    }

    vector<int64_t> flat_hub_counts(hub_idmap.size(), 0);

    for (int64_t i = 0; i < mysize; ++i)
    {
        int64_t hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            size_t loc = hub_idmap.find(hub_id)->second;
            flat_hub_counts[loc]++;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, flat_hub_counts.data(), static_cast<int>(flat_hub_counts.size()), MPI_INT64_T, MPI_SUM, comm);

    unordered_map<int64_t, int64_t> hub_counts;

    for (auto it = hub_idmap.begin(); it != hub_idmap.end(); ++it)
    {
        int64_t hub_id = it->first;
        size_t loc = it->second;
        hub_counts.insert({hub_id, flat_hub_counts[loc]});
    }

    return hub_counts;
}

void DistCoverTree::build_tree(bool verbose)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    set_times_to_zero();
    initialize_root_hub(verbose);

    double load_imbalance = numeric_limits<double>::max();
    unordered_map<int64_t, int> hub_assignments;

    while (!hub_chains.empty() && load_imbalance > 1.25)
    {
        niters++;
        compute_farthest_hub_pts(verbose);
        update_hub_chains(verbose);
        process_leaf_chains(verbose);
        process_split_chains(verbose);
        update_dists_and_pointers(verbose);
        tie(load_imbalance, hub_assignments) = compute_hub_assignments(verbose);
    }

    assert(!hub_chains.empty());
    collect_replicate_points(verbose);

    build_local_trees(hub_assignments, verbose);
}

void DistCoverTree::print_timing_results() const
{
    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        fprintf(stderr, "[tottime=%.4f,percent=%.2f] :: (initialize_root_hub)\n", initialize_root_hub_time, 100.0*(initialize_root_hub_time/overall_time));
        fprintf(stderr, "[tottime=%.4f,percent=%.2f] :: (compute_farthest_hub_pts)\n", compute_farthest_hub_pts_time, 100.0*(compute_farthest_hub_pts_time/overall_time));
        fprintf(stderr, "[tottime=%.4f,percent=%.2f] :: (update_hub_chains)\n", update_hub_chains_time, 100.0*(update_hub_chains_time/overall_time));
        fprintf(stderr, "[tottime=%.4f,percent=%.2f] :: (process_leaf_chains)\n", process_leaf_chains_time, 100.0*(process_leaf_chains_time/overall_time));
        fprintf(stderr, "[tottime=%.4f,percent=%.2f] :: (process_split_chains)\n", process_split_chains_time, 100.0*(process_split_chains_time/overall_time));
        fprintf(stderr, "[tottime=%.4f,percent=%.2f] :: (update_dists_and_pointers)\n", update_dists_and_pointers_time, 100.0*(update_dists_and_pointers_time/overall_time));
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
    nlevels = max(vertex_level+1, nlevels);

    return vertex_id;
}

double DistCoverTree::vertex_ball_radius(int64_t vertex_id) const
{
    return cover_map.find(vertex_id)->second;
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

void DistCoverTree::initialize_root_hub(bool verbose)
{
    MPITimer timer(comm, 0);
    timer.start_timer();

    my_dists.resize(mysize);
    my_hub_vtx_ids.resize(mysize);
    my_hub_pt_ids.resize(mysize);

    Point root_pt = mypoints.front();
    int64_t root_id = add_vertex(0, -1);
    hub_chains.insert({root_id, {pt[root_id]}});

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);
    MPI_Bcast(&root_pt, 1, MPI_POINT, 0, comm);
    MPI_Type_free(&MPI_POINT);

    int64_t argmax = -1;

    for (int64_t i = 0; i < mysize; ++i)
    {
        my_dists[i] = root_pt.distance(mypoints[i]);
        my_hub_vtx_ids[i] = root_id;
        my_hub_pt_ids[i] = pt[root_id];

        if (my_dists[i] > max_radius)
        {
            max_radius = my_dists[i];
            argmax = myoffset + i;
        }
    }

    MPI_Datatype MPI_ARGMAX_PAIR;
    MPI_Op MPI_ARGMAX;

    ArgmaxPair::create_mpi_dtype(&MPI_ARGMAX_PAIR);
    ArgmaxPair::create_mpi_op(&MPI_ARGMAX);

    ArgmaxPair argpair(argmax, max_radius);

    MPI_Allreduce(MPI_IN_PLACE, &argpair, 1, MPI_ARGMAX_PAIR, MPI_ARGMAX, comm);

    max_radius = argpair.value;
    argmax = argpair.index;

    timer.stop_timer();

    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        initialize_root_hub_time += timer.get_max_time();
        overall_time += timer.get_max_time();
        if (verbose) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f,itr=%lld] :: (initialize_root_hub) [argmax=%lld,max_radius=%.4f]\n", timer.get_max_time(), timer.get_avg_time(), niters, argmax, max_radius);
    }
}



void DistCoverTree::compute_farthest_hub_pts(bool verbose)
{
    /*
     * Go through all active hubs and find the point farthest from
     * the hub's current chain.
     */

    MPITimer timer(comm, 0);
    timer.start_timer();

    unordered_map<int64_t, pair<int64_t, double>> my_argmaxes; // maps hub id to (point id, distance) pair
    transform(hub_chains.begin(), hub_chains.end(), inserter(my_argmaxes, my_argmaxes.end()),
            [](const auto& chain) { return make_pair(chain.first, make_pair(-1, -1.0)); });

    // go through local points
    for (int64_t i = 0; i < mysize; ++i)
    {
        int64_t hub_id = my_hub_vtx_ids[i];

        // if point is in active hub
        if (hub_id >= 0)
        {
            auto& it = my_argmaxes.find(hub_id)->second;

            // update hub's farthest point if necessary
            if (my_dists[i] > it.second)
            {
                it.first = i;
                it.second = my_dists[i];
            }
        }
    }

    /*
     * Now each processor knows which point within its local points
     * is farthest within each hub. To determine the global farthest
     * point for each hub we do an argmax reduction across all processors
     */

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

    timer.stop_timer();

    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        compute_farthest_hub_pts_time += timer.get_max_time();
        overall_time += timer.get_max_time();
        if (verbose) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f,itr=%lld] :: (proc_chains) [hub_chains=%lu,levels=%lld]\n", timer.get_max_time(), timer.get_avg_time(), niters, hub_chains.size(), nlevels);
    }
}

void DistCoverTree::update_hub_chains(bool verbose)
{
    /*
     * Go through each hub and, based on the computed farthest point,
     * determine whether the farthest point indicates:
     *
     *     (i) that the hub is (a) a singleton or (b) a set of duplicate points -> leaf chain
     *    (ii) that the hub chain should be partitioned into new hubs -> split chain
     *   (iii) that the hub chain is incomplete -> extend chain
     *
     * Because the hubs and their farthest points are global, the routine
     * below is run identically (and redundantly) on all processors with
     * no synchronization
     */

    MPITimer timer(comm, 0);
    timer.start_timer();

    int64_t hub_id;
    pair<int64_t, double> farthest_pt;
    split_chains.clear(), leaf_chains.clear();
    int64_t extended = 0;

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
            extended++;
        }
    }

    timer.stop_timer();

    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        update_hub_chains_time += timer.get_max_time();
        overall_time += timer.get_max_time();
        if (verbose) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f,itr=%lld] :: (update_chains) [hleaves=%lu,splits=%lu,exts=%lld]\n", timer.get_max_time(), timer.get_avg_time(), niters, leaf_chains.size(), split_chains.size(), extended);
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

void DistCoverTree::process_leaf_chains(bool verbose)
{
    /*
     * Remove all leaf hubs and add associated vertices. Each
     * processor determines which of its own points are in a
     * leaf hub and accordingly removes it from the active set of points.
     *
     * Each processor batches the new added vertices it is contributing,
     * and at the end we allgather the batched vertices each processor
     * contributed and they all add them to their local copy of the tree.
     */

    MPITimer timer(comm, 0);
    timer.start_timer();

    int64_t mynlpts = 0, nlpts;

    if (!leaf_chains.empty())
    {
        for (int64_t i = 0; i < mysize; ++i)
        {
            int64_t hub_id = my_hub_vtx_ids[i];

            if (leaf_chains.find(hub_id) != leaf_chains.end())
            {
                mynlpts++;
                batch_new_vertex(i + myoffset, hub_id);
                my_hub_vtx_ids[i] = my_hub_pt_ids[i] = -1;
                my_dists[i] = 0;
            }
        }
    }

    add_batched_vertices();
    timer.stop_timer();

    if (verbose) MPI_Reduce(&mynlpts, &nlpts, 1, MPI_INT64_T, MPI_SUM, 0, comm);

    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        process_leaf_chains_time += timer.get_max_time();
        overall_time += timer.get_max_time();

        if (verbose) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f,itr=%lld] :: (proc_leaves) [leaf_pts=%lld]\n", timer.get_max_time(), timer.get_avg_time(), niters, nlpts);
    }
}

void DistCoverTree::process_split_chains(bool verbose)
{
    /*
     * Partition split hubs into new hubs, one for each point in the split chain.
     * All the points in the split hub are assigned to one of the new hubs,
     * and vertices for each new hub are added. The original split hub is deleted.
     *
     * Because the split chains are globally accessible, all processors can add
     * their new split vertices redundantly in parallel with no synchronization.
     */

    MPITimer timer(comm, 0);
    timer.start_timer();

    int64_t mynsplts = 0, nsplts;

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
            {
                mynsplts++;
                my_hub_vtx_ids[i] = it->second;
            }
        }
    }

    timer.stop_timer();

    if (verbose) MPI_Reduce(&mynsplts, &nsplts, 1, MPI_INT64_T, MPI_SUM, 0, comm);

    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        process_split_chains_time += timer.get_max_time();
        overall_time += timer.get_max_time();
        if (verbose) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f,itr=%lld] :: (proc_splits) [split_pts=%lld]\n", timer.get_max_time(), timer.get_avg_time(), niters, nsplts);
    }
}

void DistCoverTree::update_dists_and_pointers(bool verbose)
{
    /*
     * Now that hubs have been split and/or extended, go through
     * all the points and update their hub pointers and distances.
     */

    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm, 0);
    timer.start_timer();

    vector<int64_t> my_last_chain_pt_ids;
    vector<Point> my_last_chain_pts;

    for (auto it = hub_chains.begin(); it != hub_chains.end(); ++it)
    {
        int64_t last_chain_pt_id = it->second.back();

        if (myoffset <= last_chain_pt_id && last_chain_pt_id < myoffset + mysize)
        {
            my_last_chain_pt_ids.push_back(last_chain_pt_id);
            my_last_chain_pts.push_back(mypoints[last_chain_pt_id - myoffset]);
        }
    }

    vector<int> recvcounts(nprocs), displs(nprocs);
    recvcounts[myrank] = my_last_chain_pts.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    displs.front() = 0;
    partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    vector<Point> last_chain_pts(recvcounts.back() + displs.back());
    vector<int64_t> last_chain_pt_ids(recvcounts.back() + displs.back());

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);

    MPI_Allgatherv(my_last_chain_pts.data(), recvcounts[myrank], MPI_POINT, last_chain_pts.data(), recvcounts.data(), displs.data(), MPI_POINT, comm);
    MPI_Allgatherv(my_last_chain_pt_ids.data(), recvcounts[myrank], MPI_INT64_T, last_chain_pt_ids.data(), recvcounts.data(), displs.data(), MPI_INT64_T, comm);

    MPI_Type_free(&MPI_POINT);

    unordered_map<int64_t, Point> last_chain_pt_map;

    for (int64_t i = 0; i < last_chain_pts.size(); ++i)
    {
        last_chain_pt_map.insert({last_chain_pt_ids[i], last_chain_pts[i]});
    }

    for (int64_t i = 0; i < mysize; ++i)
    {
        int64_t hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            int64_t last_chain_pt_id = hub_chains.find(hub_id)->second.back();
            double lastdist = my_dists[i];
            double curdist = mypoints[i].distance(last_chain_pt_map.find(last_chain_pt_id)->second);

            if (curdist <= lastdist)
            {
                my_dists[i] = curdist;
                my_hub_pt_ids[i] = last_chain_pt_id;
            }
        }
    }

    timer.stop_timer();

    if (!myrank)
    {
        update_dists_and_pointers_time += timer.get_max_time();
        overall_time += timer.get_max_time();

        if (verbose) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f,itr=%lld] :: (updates)\n", timer.get_max_time(), timer.get_avg_time(), niters);
    }
}

vector<vector<int64_t>> DistCoverTree::build_epsilon_graph(double radius) const
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    vector<int> recvcounts(nprocs);
    recvcounts[myrank] = mysize;
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    vector<int> displs(nprocs);
    displs.front() = 0;
    partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
    assert(recvcounts.back() + displs.back() == totsize);

    vector<Point> allpoints(totsize);

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);
    MPI_Allgatherv(mypoints.data(), recvcounts[myrank], MPI_POINT, allpoints.data(), recvcounts.data(), displs.data(), MPI_POINT, comm);
    MPI_Type_free(&MPI_POINT);

    unordered_set<int64_t> idset;
    vector<int64_t> stack;
    vector<vector<int64_t>> my_edges(mysize);

    for (int64_t i = 0; i < mysize; ++i)
    {
        Point query = mypoints[i];
        idset.clear();
        stack.assign({0});

        while (!stack.empty())
        {
            int64_t u = stack.back(); stack.pop_back();

            if (query.distance(allpoints[pt[u]]) <= radius)
                idset.insert(pt[u]);

            double compare_radius = radius + max_radius * (vertex_ball_radius(u) / base);

            for (int64_t v : children[u])
                if (query.distance(allpoints[pt[v]]) <= compare_radius)
                    stack.push_back(v);
        }

        my_edges[i].reserve(idset.size());
        my_edges[i].assign(idset.begin(), idset.end());
    }

    return my_edges;
}

pair<double, unordered_map<int64_t, int>> DistCoverTree::compute_hub_assignments(bool verbose) const
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm, 0);
    timer.start_timer();

    unordered_map<int64_t, int> hub_assignments;
    vector<int64_t> workloads(nprocs, 0);
    unordered_map<int64_t, int64_t> hub_counts = get_hub_counts();

    for (auto it = hub_counts.begin(); it != hub_counts.end(); ++it)
    {
        int smallest_rank = distance(workloads.begin(), min_element(workloads.begin(), workloads.end()));
        hub_assignments.insert({it->first, smallest_rank});
        workloads[smallest_rank] += it->second;
    }

    int64_t maxcount = *max_element(workloads.begin(), workloads.end());
    int64_t totcount = accumulate(workloads.begin(), workloads.end(), 0, plus<int64_t>());
    double load_imbalance = ((nprocs + 0.0) * maxcount) / (totcount + 0.0);

    timer.stop_timer();

    if (!myrank && verbose)
    {
        fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f,itr=%lld] :: (compute_hub_assignments) :: [load_imbalance=%.4f]\n", timer.get_max_time(), timer.get_avg_time(), niters, load_imbalance);
    }

    return make_pair(load_imbalance, hub_assignments);
}

void DistCoverTree::collect_replicate_points(bool verbose)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm, 0);
    timer.start_timer();

    unordered_set<int64_t> pts(pt.begin(), pt.end());

    vector<int64_t> my_pt_ids, pt_ids;
    vector<Point> my_pt_buf, pt_buf;

    for (int64_t id : pts)
    {
        if (myoffset <= id && id < myoffset + mysize)
        {
            my_pt_ids.push_back(id);
            my_pt_buf.push_back(mypoints[id - myoffset]);
        }
    }

    vector<int> recvcounts(nprocs), displs(nprocs);
    recvcounts[myrank] = my_pt_ids.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    displs.front() = 0;
    partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    pt_ids.resize(recvcounts.back() + displs.back());
    pt_buf.resize(recvcounts.back() + displs.back());

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);

    MPI_Allgatherv(my_pt_buf.data(), recvcounts[myrank], MPI_POINT, pt_buf.data(), recvcounts.data(), displs.data(), MPI_POINT, comm);
    MPI_Allgatherv(my_pt_ids.data(), recvcounts[myrank], MPI_INT64_T, pt_ids.data(), recvcounts.data(), displs.data(), MPI_INT64_T, comm);

    MPI_Type_free(&MPI_POINT);

    for (int64_t i = 0; i < pt_buf.size(); ++i)
    {
        repoints.insert({pt_ids[i], pt_buf[i]});
    }

    timer.stop_timer();

    if (!myrank && verbose)
    {
        fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (collect_replicate_points) [num_points=%lu]\n", timer.get_max_time(), timer.get_avg_time(), repoints.size());
    }
}

void DistCoverTree::build_local_trees(const unordered_map<int64_t, int>& hub_assignments, bool verbose)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm, 0);
    timer.start_timer();

    vector<vector<int64_t>> hub_ids_sendbufs(nprocs);
    vector<vector<int64_t>> pt_ids_sendbufs(nprocs);
    vector<vector<Point>> pt_sendbufs(nprocs);
    vector<int> sendcounts(nprocs, 0), sdispls(nprocs);

    for (int64_t i = 0; i < mysize; ++i)
    {
        int64_t hub_id = my_hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            if (pt[hub_id] != i + myoffset)
            {
                int rank = hub_assignments.find(hub_id)->second;
                hub_ids_sendbufs[rank].push_back(hub_id);
                pt_ids_sendbufs[rank].push_back(myoffset + i);
                pt_sendbufs[rank].push_back(mypoints[i]);
                sendcounts[rank]++;
            }
        }
    }

    sdispls.front() = 0;
    partial_sum(sendcounts.begin(), sendcounts.end()-1, sdispls.begin()+1);

    vector<int64_t> hub_ids_sendbuf;
    vector<int64_t> pt_ids_sendbuf;
    vector<Point> pt_sendbuf;

    for (int i = 0; i < nprocs; ++i)
    {
        hub_ids_sendbuf.insert(hub_ids_sendbuf.end(), hub_ids_sendbufs[i].begin(), hub_ids_sendbufs[i].end());
        pt_ids_sendbuf.insert(pt_ids_sendbuf.end(), pt_ids_sendbufs[i].begin(), pt_ids_sendbufs[i].end());
        pt_sendbuf.insert(pt_sendbuf.end(), pt_sendbufs[i].begin(), pt_sendbufs[i].end());
    }

    vector<int> recvcounts(nprocs), rdispls(nprocs);

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    rdispls.front() = 0;
    partial_sum(recvcounts.begin(), recvcounts.end()-1, rdispls.begin()+1);

    vector<Point> pt_recvbuf(recvcounts.back() + rdispls.back());
    vector<int64_t> hub_ids_recvbuf(recvcounts.back() + rdispls.back());
    vector<int64_t> pt_ids_recvbuf(recvcounts.back() + rdispls.back());

    MPI_Datatype MPI_POINT;
    Point::create_mpi_dtype(&MPI_POINT);

    MPI_Alltoallv(pt_sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_POINT,
                  pt_recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_POINT, comm);

    MPI_Alltoallv(hub_ids_sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT64_T,
                  hub_ids_recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_INT64_T, comm);

    MPI_Alltoallv(pt_ids_sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT64_T,
                  pt_ids_recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_INT64_T, comm);

    MPI_Type_free(&MPI_POINT);

    local_ptid_map.clear(), local_trees.clear();
    unordered_map<int64_t, vector<Point>> local_pt_map;
    unordered_set<int64_t> hub_idset(hub_ids_recvbuf.begin(), hub_ids_recvbuf.end());

    for (int64_t hub_id : hub_idset)
    {
        local_ptid_map.insert({hub_id, {pt[hub_id]}});
        local_pt_map.insert({hub_id, {repoints.find(pt[hub_id])->second}});
    }

    for (int64_t i = 0; i < pt_recvbuf.size(); ++i)
    {
        Point p = pt_recvbuf[i];
        int64_t hub_id = hub_ids_recvbuf[i];
        int64_t pt_id = pt_ids_recvbuf[i];

        vector<int64_t>& ptids = local_ptid_map.find(hub_id)->second;
        vector<Point>& pts = local_pt_map.find(hub_id)->second;

        ptids.push_back(pt_id);
        pts.push_back(p);
    }

    for (auto it = local_ptid_map.begin(); it != local_ptid_map.end(); ++it)
    {
        local_trees.insert({it->first, CoverTree(local_pt_map.find(it->first)->second, base)});
    }

    for (auto it = local_trees.begin(); it != local_trees.end(); ++it)
    {
        it->second.build_tree();
    }

    timer.stop_timer();

    vector<int64_t> num_trees;
    if (!myrank) num_trees.resize(nprocs);
    int64_t my_num_trees = local_trees.size();
    MPI_Gather(&my_num_trees, 1, MPI_INT64_T, num_trees.data(), 1, MPI_INT64_T, 0, comm);

    if (verbose && !myrank)
    {
        fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (build_local_trees)\n", timer.get_max_time(), timer.get_avg_time());

        stringstream ss;
        for (int i = 0; i < nprocs; ++i) ss << "[myrank=" << i << ",num_trees=" << num_trees[i] << "]\n";
        string s = ss.str();
        fprintf(stderr, "%s", s.c_str());
        fflush(stderr);
    }
}
