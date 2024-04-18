#include "CoverTree.h"
#include "Point.h"
#include <limits>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <iomanip>
#include <cassert>
#include <stdio.h>

CoverTree::CoverTree(const vector<Point>& points, double base)
    : max_radius(-1), base(base), points(points)
{
    build_tree_point_loop();
    //build_tree_hub_loop();
}

vector<int64_t> CoverTree::radii_query(const Point& query, double radius) const
{
    unordered_set<int64_t> idset; // point ids of points within a distance of @radius of @query
    vector<int64_t> stack = {0}; // start exploring tree from root (point id 0)

    while (stack.size() != 0)
    {
        // get next point @u
        int64_t u = stack.back(); stack.pop_back();
        vector<int64_t> mychildren = children[u];

        // go through children of @u
        for (int64_t v : mychildren)
        {
            // if child of @u is within a distance of @radius and within a distance
            // the covering distance of @u then we need to explore its children
            if (query.distance(get_vertex_point(v)) <= radius + max_radius*vertex_ball_radius(v))
                stack.push_back(v);

            // add @radius-neighbor of @u to set
            if (query.distance(get_vertex_point(v)) <= radius)
                idset.insert(pt[v]);
        }
    }

    return vector<int64_t>(idset.begin(), idset.end());
}

vector<vector<int64_t>> CoverTree::build_epsilon_graph(double radius) const
{
    // graph adjacency list (initialized empty)
    vector<vector<int64_t>> graph(num_points());

    for (int64_t u = 0; u < num_points(); ++u)
    {
        // compute @radius-neighbors of @u
        graph[u] = radii_query(get_point(u), radius);
    }

    return graph;
}

unordered_map<int64_t, int64_t> find_argmaxes(const vector<double>& dists, const vector<int64_t>& hub_vtx_ids, const unordered_map<int64_t, vector<int64_t>>& hub_chains)
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

    unordered_map<int64_t, int64_t> argmax_ids;
    transform(argmaxes.begin(), argmaxes.end(), inserter(argmax_ids, argmax_ids.end()),
              [](auto pair) { return make_pair(pair.first, pair.second.first); });

    return argmax_ids;
}

/*
 * CoverTree::set_max_radius() determines the maximum distance between the
 * root point (the first point in the data by convention) and all other
 * points. This is used to normalize distance calculations so that the
 * maximum normalized distance between the root vertex and all other vertices
 * in the tree is 1.
 */
void CoverTree::set_max_radius()
{
    // to normalize distances, find distance of farthest point from
    // root (vertex index 0)
    for (const Point& point : points)
    {
        max_radius = max(max_radius, points.front().distance(point));
    }
}

void CoverTree::build_tree_hub_loop()
{
    set_max_radius();

    // list of tasks, each of which is identified by a (1) unexpanded vertex (newly constructed)
    // and (2) list of its descendants.
    list<tuple<int64_t, vector<int64_t>>> tasks;

    // list of descendant points for vertex we are currently "expanding"
    vector<int64_t> descendants(num_points());
    iota(descendants.begin(), descendants.end(), 0); // starts with all of them

    // add the root task with all points as descendants
    int64_t parent_id = add_vertex(descendants.front(), -1);
    tasks.emplace_back(parent_id, descendants);

    // for the current task (expanding a vertex with a given vertex id)
    // we partition the set of descendants of that vertex into "Voronoi cells",
    // with the first point id in each partition being a point in the "Voronoi
    // diagram" and every other point in a given partition will be a tree
    // descendant of the first point
    vector<vector<int64_t>> child_partitions;

    // while there are still unexpanded vertices in tree
    while (tasks.size() != 0)
    {
        // get next vertex to expand and its associated descendant points
        tie(parent_id, descendants) = tasks.front(); tasks.pop_front();

        // compute a partitioning of @descendants ids so that the first
        // id of each partition will be a direct child of the vertex @parent_id
        // and so that the set of all descendant points that are closest to
        // the direct child comprise the remaining points in the partition
        child_partitions = compute_child_partitions(parent_id, descendants);

        // add all direct children as vertices in the graph and create new tasks
        // for each of them with their associated descendants partition.
        for (int64_t i = 0; i < child_partitions.size(); ++i)
        {
            int64_t vertex_id = add_vertex(child_partitions[i].front(), parent_id);
            tasks.emplace_back(vertex_id, child_partitions[i]);
        }
    }
}

void CoverTree::build_tree_point_loop()
{
    vector<double> dists(num_points());
    vector<int64_t> hub_vtx_ids(num_points());
    vector<int64_t> hub_pt_ids(num_points());

    int64_t root_id = add_vertex(0, -1);

    for (int64_t i = 0; i < num_points(); ++i)
    {
        dists[i] = get_vertex_point(root_id).distance(get_point(i));
        hub_vtx_ids[i] = root_id;
        hub_pt_ids[i] = pt[root_id];
        max_radius = max(dists[i], max_radius);
    }

    unordered_map<int64_t, vector<int64_t>> hub_chains;
    hub_chains.insert({root_id, {pt[root_id]}});

    while (hub_chains.size() != 0)
    {
        unordered_map<int64_t, int64_t> argmaxes = find_argmaxes(dists, hub_vtx_ids, hub_chains);

        int64_t hub_id, farthest_pt_id;
        unordered_set<int64_t> stopped_chains, leaf_chains;

        for (auto it = argmaxes.begin(); it != argmaxes.end(); ++it)
        {
            tie(hub_id, farthest_pt_id) = *it;
            double farthest_dist = dists[farthest_pt_id] / max_radius;

            if (farthest_dist == 0)
            {
                for (int64_t i = 0; i < num_points(); ++i)
                    if (hub_id == hub_vtx_ids[i])
                        add_vertex(i, hub_id);

                /*
                 * TODO: The above loop is pretty inefficient here, but temporarily necessary to deal
                 * with the problem of duplicate points. If there were no duplicate points, then
                 * we could get away with the following:
                 *
                 * """
                 * const vector<int64_t>& leaves = hub_chains.find(hub_id)->second;
                 *
                 * if (leaves.size() != 1)
                 *     for (int64_t leaf_pt : leaves)
                 *         add_vertex(leaf_pt, hub_id);
                 * """
                 *
                 * However, because there may be duplicate points we need to make sure that
                 * we are not neglecting the possibility that all the hub points are identical
                 * to the hub parent vertex, in which case the just mentioned alternative
                 * solution would miss all the duplicates.
                 */

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

        if (leaf_chains.size() != 0)
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

        if (stopped_chains.size() != 0)
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
}

void CoverTree::update_cell(vector<double>& dists, vector<int64_t>& closest, int64_t farthest, const vector<int64_t>& descendants) const
{
    int64_t n = dists.size();
    assert(dists.size() == closest.size());

    // compute the distance between @descendants[@farthest] and every other
    // point in @descendants and update @dists with the minimum of whatever
    // value was computed and the previous value that was stored there; also
    // update @closest vector to point to @descendants[@farthest] for those
    // descendant points that became closer.
    for (int64_t i = 0; i < n; ++i)
    {
        double lastdist = dists[i];
        double curdist = get_point(descendants[i]).distance(get_point(descendants[farthest])) / max_radius;

        if (curdist <= lastdist)
        {
            dists[i] = curdist;
            closest[i] = descendants[farthest];
        }
    }
}

vector<vector<int64_t>> CoverTree::compute_child_partitions(int64_t parent_id, const vector<int64_t>& descendants) const
{
    // @child_partitions simply partitions @descendants into disjoint sets of
    // direct children (first point id in partition vector) and their descendants
    // (remaining points in partition vector)
    vector<vector<int64_t>> child_partitions;
    int64_t n_descendants = descendants.size();

    // first descendant point must be the nested point of parent
    assert(n_descendants >= 1 && pt[parent_id] == descendants.front());

    // if no descendants (besides nested child) then @parent_id
    // is a leaf vertex and we don't compute any more of its descendants
    if (n_descendants == 1) return child_partitions;

    // maintain distances between growing set of computed direct children
    // and the set of all @descendants, in addition to pointers telling us
    // each of @descendants closest direct children.
    vector<double> dists(n_descendants, numeric_limits<double>::max());
    vector<int64_t> closest(n_descendants, descendants.front());

    // initialize distances and pointers of all @descendants relative to
    // the parent point (same as nested point)
    update_cell(dists, closest, 0, descendants);

    // if distance between parent point and all of its descendants is 0,
    // then it means that all of these points are duplicates; therefore
    // we just add them all as direct children leaves of parent and no
    // more descendants need to be computed of @parent_id
    if (all_of(dists.begin(), dists.end(), [](double d) { return d == 0; }))
    {
        for (int64_t duplicate_child_pt : descendants)
        {
            child_partitions.push_back({duplicate_child_pt});
        }
        return child_partitions;
    }

    // initialize nested child partition hub
    child_partitions.push_back({pt[parent_id]});

    // cover radius of parent for determining direct children
    double parent_ball_radius = vertex_ball_radius(parent_id);

    // loop through descendants
    for (int64_t k = 0; k < n_descendants; ++k)
    {
        // find descendant farthest from current set of partition hub points
        int64_t farthest = distance(dists.begin(), max_element(dists.begin(), dists.end()));

        // if we violate separation condition here then stop adding
        // new direct children
        if (dists[farthest] <= (parent_ball_radius / base))
            break;

        // add new partition and associated hub point
        child_partitions.push_back({descendants[farthest]});

        // update distances and pointers
        update_cell(dists, closest, farthest, descendants);
    }

    // go through each partition and use pointers to determine
    // which non-child descendants are closest to partition hub point
    // and add them to partition; don't add partition hub point itself
    // since we already added it (because it is assumed to always be first
    // in partition)
    for (int64_t i = 0; i < child_partitions.size(); ++i)
    {
        int64_t hub_point = child_partitions[i].front();

        for (int64_t j = 0; j < n_descendants; ++j)
            if (closest[j] == hub_point && descendants[j] != hub_point)
                child_partitions[i].push_back(descendants[j]);
    }

    return child_partitions;
}

/*
 * CoverTree::add_vertex(@point_id, @parent_id) adds a vertex to
 * the tree whose associated point has point id @point_id and whose
 * parent vertex has vertex id @parent_id.
 */
int64_t CoverTree::add_vertex(int64_t point_id, int64_t parent_id)
{
    int64_t vertex_level; // level of new vertex
    int64_t vertex_id = pt.size(); // vertex id of new vertex

    pt.push_back(point_id); // pt[vertex_id] = point_id
    children.emplace_back(); // new vertex initialized with no children

    if (parent_id >= 0) // @vertex_id is not the root vertex
    {
        vertex_level = level[parent_id] + 1;
        children[parent_id].push_back(vertex_id);
    }
    else
    {
        vertex_level = 0;
    }

    level.push_back(vertex_level);

    return vertex_id;
}

/*
 * CoverTree::vertex_ball_radius(@vertex_id) returns the radius
 * of the coverage for the children of @vertex_id.
 *
 * For a vertex on level k, its children must all be within a distance
 * of 1 / 2^k.
 */
double CoverTree::vertex_ball_radius(int64_t vertex_id) const
{
    return pow(base, -1. * level[vertex_id]);
}

void CoverTree::write_gml(const char *filename) const
{
    FILE *f = fopen(filename, "w");

    fprintf(f, "graph\n[\n");

    for (int64_t id = 0; id < num_vertices(); ++id)
    {
        const Point& point = get_vertex_point(id);
        fprintf(f, "\tnode\n\t[\n\t\tid %lld\n\t\tpt %lld\n\t\tlevel %lld\n\t\tcover %.3f\n\t]\n",
                id, pt[id], level[id], vertex_ball_radius(id));
    }

    vector<int64_t> stack = {0};

    while (stack.size() != 0)
    {
        int64_t vtx = stack.back(); stack.pop_back();
        vector<int64_t> mychildren = children[vtx];

        const Point& source = get_vertex_point(vtx);

        for (int64_t child : mychildren)
        {
            const Point& target = get_vertex_point(child);
            fprintf(f, "\tedge\n\t[\n\t\tsource %lld\n\t\ttarget %lld\n\t\tdistance %.3f\n\t]\n",
                    vtx, child, source.distance(target) / max_radius);
            stack.push_back(child);
        }
    }

    fprintf(f, "]\n");

    fclose(f);
}
