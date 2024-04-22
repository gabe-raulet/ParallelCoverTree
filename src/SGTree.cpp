#include "SGTree.h"
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

SGTree::SGTree(const vector<Point>& points, double base) : max_radius(-1), base(base), points(points) {}

int64_t SGTree::num_points() const { return points.size(); }
int64_t SGTree::num_vertices() const { return pt.size(); }
int64_t SGTree::num_levels() const { return *max_element(level.begin(), level.end()) + 1; }

Point SGTree::get_point(int64_t point_id) const { return points[point_id]; }
Point SGTree::get_vertex_point(int64_t vertex_id) const { return points[pt[vertex_id]]; }

int64_t SGTree::add_vertex(int64_t point_id, int64_t parent_id)
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
    return vertex_id;
}

double SGTree::vertex_ball_radius(int64_t vertex_id) const
{
    return pow(base, -1. * level[vertex_id]);
}

void SGTree::initialize_root_hub()
{
    auto t1 = mytimer::clock::now();

    dists.resize(num_points());
    hub_vtx_ids.resize(num_points());
    hub_pt_ids.resize(num_points());

    int64_t root_id = add_vertex(0, -1);

    for (int64_t i = 0; i < num_points(); ++i)
    {
        dists[i] = get_vertex_point(root_id).distance(get_point(i));
        hub_vtx_ids[i] = root_id;
        hub_pt_ids[i] = pt[root_id];
        max_radius = max(dists[i], max_radius);
    }

    hub_chains.insert({root_id, {pt[root_id]}});

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    initialize_root_hub_times.push_back(t);
}

void SGTree::compute_farthest_hub_pts()
{
    auto t1 = mytimer::clock::now();

    unordered_map<int64_t, pair<int64_t, double>> argmaxes;
    transform(hub_chains.begin(), hub_chains.end(), inserter(argmaxes, argmaxes.end()),
              [](auto pair) { return make_pair(pair.first, make_pair(-1, -1.0)); });

    for (int64_t i = 0; i < num_points(); ++i)
    {
        int64_t hub_id = hub_vtx_ids[i];

        if (hub_id >= 0)
        {
            auto& it = argmaxes.find(hub_id)->second;

            if (dists[i] > it.second)
            {
                it.first = i;
                it.second = dists[i];
            }
        }
    }

    farthest_hub_pts.clear();
    transform(argmaxes.begin(), argmaxes.end(), inserter(farthest_hub_pts, farthest_hub_pts.end()),
              [](auto pair) { return make_pair(pair.first, pair.second.first); });

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    compute_farthest_hub_pts_times.push_back(t);
}

void SGTree::update_hub_chains()
{
    auto t1 = mytimer::clock::now();

    int64_t hub_id, farthest_pt_id;
    split_chains.clear(), leaf_chains.clear();

    for (auto it = farthest_hub_pts.begin(); it != farthest_hub_pts.end(); ++it)
    {
        tie(hub_id, farthest_pt_id) = *it;
        double farthest_dist = dists[farthest_pt_id] / max_radius;

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

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    update_hub_chains_times.push_back(t);
}

void SGTree::process_leaf_chains()
{
    auto t1 = mytimer::clock::now();

    if (!leaf_chains.empty())
    {
        for (int64_t i = 0; i < num_points(); ++i)
        {
            int64_t hub_id = hub_vtx_ids[i];

            if (leaf_chains.find(hub_id) != leaf_chains.end())
            {
                add_vertex(i, hub_id);
                hub_vtx_ids[i] = hub_pt_ids[i] = -1;
                dists[i] = 0;
            }
        }
    }

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    process_leaf_chains_times.push_back(t);
}

void SGTree::process_split_chains()
{
    auto t1 = mytimer::clock::now();

    if (!split_chains.empty())
    {
        unordered_map<int64_t, int64_t> hub_pt_id_updates;
        vector<pair<int64_t, int64_t>> new_hubs;

        for (int64_t hub_id : split_chains)
        {
            const vector<int64_t>& chain = hub_chains.find(hub_id)->second;

            for (int64_t pt_id : chain)
            {
                new_hubs.emplace_back(pt_id, hub_id);
            }

            hub_chains.erase(hub_id);
        }

        for (const auto& hub_pair : new_hubs)
        {
            int64_t vtx_id = add_vertex(hub_pair.first, hub_pair.second);
            hub_chains.insert({vtx_id, {hub_pair.first}});
            hub_pt_id_updates.insert({hub_pair.first, vtx_id});
        }

        for (int64_t i = 0; i < num_points(); ++i)
        {
            int64_t closest_pt_id = hub_pt_ids[i];
            auto it = hub_pt_id_updates.find(closest_pt_id);

            if (it != hub_pt_id_updates.end())
                hub_vtx_ids[i] = it->second;
        }
    }

    auto t2 = mytimer::clock::now();
    double t = mytimer::duration(t2-t1).count();
    process_split_chains_times.push_back(t);
}

void SGTree::update_dists_and_pointers()
{
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
}

void SGTree::build_tree()
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

void SGTree::print_timing_results() const
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
