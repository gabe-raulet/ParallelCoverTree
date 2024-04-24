#include "CoverTree.h"
#include "Point.h"
#include "MyTimer.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <limits>
#include <unordered_set>
#include <iomanip>
#include <cassert>
#include <stdio.h>

CoverTree::CoverTree(const vector<Point>& points, double base) : max_radius(-1), base(base), points(points), nlevels(0), niters(0) {}

int64_t CoverTree::num_points() const { return points.size(); }
int64_t CoverTree::num_vertices() const { return pt.size(); }
int64_t CoverTree::num_levels() const { return *max_element(level.begin(), level.end()) + 1; }

Point CoverTree::get_point(int64_t point_id) const { return points[point_id]; }
Point CoverTree::get_vertex_point(int64_t vertex_id) const { return get_point(pt[vertex_id]); }

int64_t CoverTree::add_vertex(int64_t point_id, int64_t parent_id)
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

    nlevels = max(vertex_level+1, nlevels);

    level.push_back(vertex_level);
    return vertex_id;
}

double CoverTree::vertex_ball_radius(int64_t vertex_id) const
{
    return pow(base, -1. * level[vertex_id]);
}

void CoverTree::initialize_root_hub(bool verbose)
{
    auto t1 = mytimer::clock::now();

    dists.resize(num_points());
    hub_vtx_ids.resize(num_points());
    hub_pt_ids.resize(num_points());

    int64_t root_id = add_vertex(0, -1);
    int64_t argmax = -1;

    for (int64_t i = 0; i < num_points(); ++i)
    {
        dists[i] = get_vertex_point(root_id).distance(get_point(i));
        hub_vtx_ids[i] = root_id;
        hub_pt_ids[i] = pt[root_id];

        if (dists[i] > max_radius)
        {
            max_radius = dists[i];
            argmax = i;
        }
    }

    hub_chains.insert({root_id, {pt[root_id]}});

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    initialize_root_hub_times.push_back(t);

    if (verbose) fprintf(stderr, "[time=%.4f,itr=%lld] :: (init_root_hub) [argmax=%lld,max_radius=%.4f]\n", t, niters, argmax, max_radius);
}

void CoverTree::compute_farthest_hub_pts(bool verbose)
{
    /*
     * Go through all active hubs and find the point farthest from
     * the hub's current chain.
     */

    auto t1 = mytimer::clock::now();

    farthest_hub_pts.clear();
    transform(hub_chains.begin(), hub_chains.end(), inserter(farthest_hub_pts, farthest_hub_pts.end()),
              [](auto pair) { return make_pair(pair.first, make_pair(-1, -1.0)); });

    // go through points
    for (int64_t i = 0; i < num_points(); ++i)
    {
        int64_t hub_id = hub_vtx_ids[i];

        // if point is in an active hub
        if (hub_id >= 0)
        {
            auto& it = farthest_hub_pts.find(hub_id)->second;

            // update hub's farthest point if necessary
            if (dists[i] > it.second)
            {
                it.first = i;
                it.second = dists[i];
            }
        }
    }

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    compute_farthest_hub_pts_times.push_back(t);

    if (verbose) fprintf(stderr, "[time=%.4f,itr=%lld] :: (proc_chains) [hub_chains=%lu,levels=%lld]\n", t, niters, hub_chains.size(), num_levels());
}

void CoverTree::update_hub_chains(bool verbose)
{
    /*
     * Go through each hub and, based on the computed farthest point,
     * determine whether the farthest point indicates:
     *
     *     (i) that the hub is (a) a singleton or (b) a set of duplicate points -> leaf chain
     *    (ii) that the hub chain should be partitioned into new hubs -> split chain
     *   (iii) that the hub chain is incomplete -> extend chain
     */

    auto t1 = mytimer::clock::now();

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

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    update_hub_chains_times.push_back(t);

    if (verbose) fprintf(stderr, "[time=%.4f,itr=%lld] :: (update_chains) [hleaves=%lu,splits=%lu,exts=%lld]\n", t, niters, leaf_chains.size(), split_chains.size(), extended);
}

void CoverTree::process_leaf_chains(bool verbose)
{
    /*
     * Remove all leaf hubs and add associated vertices
     */

    auto t1 = mytimer::clock::now();
    int64_t nlpts = 0;

    if (!leaf_chains.empty())
    {
        for (int64_t i = 0; i < num_points(); ++i)
        {
            int64_t hub_id = hub_vtx_ids[i];

            if (leaf_chains.find(hub_id) != leaf_chains.end())
            {
                nlpts++;
                add_vertex(i, hub_id);
                hub_vtx_ids[i] = hub_pt_ids[i] = -1;
                dists[i] = 0;
            }
        }
    }

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    process_leaf_chains_times.push_back(t);

    if (verbose) fprintf(stderr, "[time=%.4f,itr=%lld] :: (proc_leaves) [leaf_pts=%lld]\n", t, niters, nlpts);
}

void CoverTree::process_split_chains(bool verbose)
{
    /*
     * Partition split hubs into new hubs, one for each point in the split chain.
     * All the points in the split hub are assigned to one of the new hubs,
     * and vertices for each new hub are added. The original split hub is deleted.
     */

    auto t1 = mytimer::clock::now();

    int64_t nsplts = 0;

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

        // add new hubs and vertices
        for (int64_t i = 0; i < new_hub_pts.size(); ++i)
        {
            int64_t vtx_id = add_vertex(new_hub_pts[i], new_hub_ids[i]);
            hub_chains.insert({vtx_id, {new_hub_pts[i]}});
            hub_pt_id_updates.insert({new_hub_pts[i], vtx_id});
        }

        // reassign points in original hub to one of the new hubs
        for (int64_t i = 0; i < num_points(); ++i)
        {
            int64_t closest_pt_id = hub_pt_ids[i];
            auto it = hub_pt_id_updates.find(closest_pt_id);

            if (it != hub_pt_id_updates.end())
            {
                nsplts++;
                hub_vtx_ids[i] = it->second;
            }
        }
    }

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    process_split_chains_times.push_back(t);

    if (verbose) fprintf(stderr, "[time=%.4f,itr=%lld] :: (proc_splits) [split_pts=%lld]\n", t, niters, nsplts);
}

void CoverTree::update_dists_and_pointers(bool verbose)
{
    /*
     * Now that hubs have been split and/or extended, go through
     * all the points and update their hub pointers and distances
     */

    auto t1 = mytimer::clock::now();

    for (int64_t i = 0; i < num_points(); ++i)
    {
        int64_t hub_id = hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            int64_t last_chain_pt_id = hub_chains.find(hub_id)->second.back();
            double lastdist = dists[i];
            double curdist = get_point(i).distance(get_point(last_chain_pt_id));

            if (curdist <= lastdist)
            {
                dists[i] = curdist;
                hub_pt_ids[i] = last_chain_pt_id;
            }
        }
    }

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    update_dists_and_pointers_times.push_back(t);

    if (verbose) fprintf(stderr, "[time=%.4f,itr=%lld] :: (updates)\n", t, niters);
}

void CoverTree::build_tree(bool verbose)
{
    initialize_root_hub(verbose);

    while (!hub_chains.empty())
    {
        niters++;
        compute_farthest_hub_pts(verbose);
        update_hub_chains(verbose);
        process_leaf_chains(verbose);
        process_split_chains(verbose);
        update_dists_and_pointers(verbose);

        if (verbose)
        {
            double tc = compute_farthest_hub_pts_times.back();
            double tu = update_hub_chains_times.back();
            double tl = process_leaf_chains_times.back();
            double ts = process_split_chains_times.back();
            double td = update_dists_and_pointers_times.back();
            double tall = tc + tu + tl + ts + td;
            fprintf(stderr, "[time=%.4f,itr=%lld] :: time splits :: [%.1f,%.1f,%.1f,%.1f,%.1f]\n", tall, niters, 100.0*(tc/tall), 100.0*(tu/tall), 100.0*(tl/tall), 100.0*(ts/tall), 100.0*(td/tall));
        }
    }
}

tuple<double, double, double> get_stats(const vector<double>& times)
{
    double tot = accumulate(times.cbegin(), times.cend(), 0.0, plus<double>());
    double avg = tot / static_cast<double>(times.size());

    vector<double> devs;
    devs.reserve(times.size());

    transform(times.cbegin(), times.cend(), back_inserter(devs), [&avg](double val) { return (val-avg) * (val-avg); });
    double var = accumulate(devs.cbegin(), devs.cend(), 0.0, plus<double>()) / static_cast<double>(devs.size());

    return {tot, avg, sqrt(var)};
}

void CoverTree::print_timing_results() const
{
    double overall = 0.0;
    overall += std::get<0>(get_stats(initialize_root_hub_times));
    overall += std::get<0>(get_stats(compute_farthest_hub_pts_times));
    overall += std::get<0>(get_stats(update_hub_chains_times));
    overall += std::get<0>(get_stats(process_leaf_chains_times));
    overall += std::get<0>(get_stats(process_split_chains_times));
    overall += std::get<0>(get_stats(update_dists_and_pointers_times));

    double tot, avg, stddev;

    tie(tot, avg, stddev) = get_stats(initialize_root_hub_times);
    fprintf(stderr, "[tottime=%.4f,avgtime=%.4f,sdtime=%.4f,percent=%.4f] :: (initialize_root_hub)\n", tot, avg, stddev, 100.0*(tot/overall));

    tie(tot, avg, stddev) = get_stats(compute_farthest_hub_pts_times);
    fprintf(stderr, "[tottime=%.4f,avgtime=%.4f,sdtime=%.4f,percent=%.4f] :: (compute_farthest_hub_pts)\n", tot, avg, stddev, 100.0*(tot/overall));

    tie(tot, avg, stddev) = get_stats(update_hub_chains_times);
    fprintf(stderr, "[tottime=%.4f,avgtime=%.4f,sdtime=%.4f,percent=%.4f] :: (update_hub_chains)\n", tot, avg, stddev, 100.0*(tot/overall));

    tie(tot, avg, stddev) = get_stats(process_leaf_chains_times);
    fprintf(stderr, "[tottime=%.4f,avgtime=%.4f,sdtime=%.4f,percent=%.4f] :: (process_leaf_chains)\n", tot, avg, stddev, 100.0*(tot/overall));

    tie(tot, avg, stddev) = get_stats(process_split_chains_times);
    fprintf(stderr, "[tottime=%.4f,avgtime=%.4f,sdtime=%.4f,percent=%.4f] :: (process_split_chains)\n", tot, avg, stddev, 100.0*(tot/overall));

    tie(tot, avg, stddev) = get_stats(update_dists_and_pointers_times);
    fprintf(stderr, "[tottime=%.4f,avgtime=%.4f,sdtime=%.4f,percent=%.4f] :: (update_dists_and_pointers)\n", tot, avg, stddev, 100.0*(tot/overall));
}

vector<int64_t> CoverTree::radii_query(const Point& query, double radius) const
{
    unordered_set<int64_t> idset;
    vector<int64_t> stack = {0};

    while (stack.size() != 0)
    {
        int64_t u = stack.back(); stack.pop_back();
        vector<int64_t> mychildren = children[u];

        for (int64_t v : mychildren)
        {
            if (query.distance(get_vertex_point(v)) <= radius + max_radius*vertex_ball_radius(v))
                stack.push_back(v);

            if (query.distance(get_vertex_point(v)) <= radius)
                idset.insert(pt[v]);
        }
    }

    return vector<int64_t>(idset.begin(), idset.end());
}

vector<vector<int64_t>> CoverTree::build_epsilon_graph(double radius) const
{
    vector<vector<int64_t>> graph(num_points());

    for (int64_t u = 0; u < num_points(); ++u)
    {
        graph[u] = radii_query(get_point(u), radius);
    }

    return graph;
}
