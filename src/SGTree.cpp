#include "SGTree.h"
#include "Point.h"
#include "MyTimer.h"
#include <iostream>
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
}

void SGTree::compute_farthest_hub_pts()
{
    unordered_map<int64_t, pair<int64_t, double>> argmaxes;
    transform(hub_chains.begin(), hub_chains.end(), inserter(argmaxes, argmaxes.end()),
              [](auto pair) { return make_pair(pair.first, make_pair(-1, -1.0)); });

    for (int64_t i = 0; i < dists.size(); ++i)
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
}

void SGTree::update_hub_chains()
{
    int64_t hub_id, farthest_pt_id;
    stopped_chains.clear(), leaf_chains.clear();

    for (auto it = farthest_hub_pts.begin(); it != farthest_hub_pts.end(); ++it)
    {
        tie(hub_id, farthest_pt_id) = *it;
        double farthest_dist = dists[farthest_pt_id] / max_radius;

        if (farthest_dist == 0)
        {
            for (int64_t i = 0; i < num_points(); ++i)
                if (hub_id == hub_vtx_ids[i])
                    add_vertex(i, hub_id);

            hub_chains.erase(hub_id);
            leaf_chains.insert(hub_id);
        }
        else if (farthest_dist <= (vertex_ball_radius(hub_id) / base))
        {
            stopped_chains.insert(hub_id);
        }
        else
        {
            hub_chains.find(hub_id)->second.push_back(farthest_pt_id);
        }
    }
}

void SGTree::process_leaf_chains()
{
    if (!leaf_chains.empty())
    {
        for (int64_t i = 0; i < num_points(); ++i)
        {
            int64_t hub_id = hub_vtx_ids[i];

            if (leaf_chains.find(hub_id) != leaf_chains.end())
            {
                hub_vtx_ids[i] = hub_pt_ids[i] = -1;
                dists[i] = 0;
            }
        }
    }
}

void SGTree::process_split_chains()
{
    if (!stopped_chains.empty())
    {
        unordered_map<int64_t, int64_t> hub_pt_id_updates;

        for (int64_t hub_id : stopped_chains)
        {
            const vector<int64_t>& chain = hub_chains.find(hub_id)->second;

            for (int64_t pt_id : chain)
            {
                int64_t vtx_id = add_vertex(pt_id, hub_id);
                hub_chains.insert({vtx_id, {pt_id}});
                hub_pt_id_updates.insert({pt_id, vtx_id});
            }

            hub_chains.erase(hub_id);
        }

        for (int64_t i = 0; i < num_points(); ++i)
        {
            int64_t closest_pt_id = hub_pt_ids[i];
            auto it = hub_pt_id_updates.find(closest_pt_id);

            if (it != hub_pt_id_updates.end())
                hub_vtx_ids[i] = it->second;
        }
    }
}

void SGTree::update_dists_and_pointers()
{
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
