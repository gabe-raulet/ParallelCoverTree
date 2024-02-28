#include "CoverTree.h"
#include <cmath>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <omp.h>

double CoverTree::point_dist(int64_t point_id1, int64_t point_id2) const
{
    return distance(points[point_id1], points[point_id2]) / max_radius;
}

double CoverTree::vertex_ball_radius(int64_t vertex_id) const
{
    return pow(2., -1. * get_vertex_level(vertex_id));
}

size_t CoverTree::get_child_ids(int64_t vertex_id, std::vector<int64_t>& child_ids) const
{
    child_ids.clear();
    const auto& ids = children[vertex_id];
    std::copy(ids.begin(), ids.end(), std::back_inserter(child_ids));
    return child_ids.size();
}

size_t CoverTree::get_level_ids(int64_t vertex_level, std::vector<int64_t>& level_ids) const
{
    level_ids.clear();
    const auto& ids = levelset[vertex_level];
    std::copy(ids.begin(), ids.end(), std::back_inserter(level_ids));
    return level_ids.size();
}

size_t CoverTree::radii_query(const Point& query, double radius, std::vector<int64_t>& ids) const
{
    std::unordered_set<int64_t> idset;
    std::vector<int64_t> stack = {0};

    int64_t u, v;
    std::vector<int64_t> mychildren;
    size_t m;

    while (stack.size() != 0)
    {
        u = stack.back(); stack.pop_back();
        m = get_child_ids(u, mychildren);

        for (size_t i = 0; i < m; ++i)
        {
            v = mychildren[i];

            if (distance(query, points[get_point_id(v)]) <= radius + max_radius*vertex_ball_radius(v))
            {
                stack.push_back(v);
            }

            if (distance(query, points[get_point_id(v)]) <= radius)
            {
                idset.insert(get_point_id(v));
            }
        }
    }

    ids.assign(idset.begin(), idset.end());
    return ids.size();
}

int64_t CoverTree::add_vertex(int64_t point_id, int64_t parent_id)
{
    int64_t vertex_level;
    int64_t vertex_id = num_vertices();

    pt.push_back(point_id);
    parent.push_back(parent_id);
    children.resize(children.size()+1);

    if (parent_id >= 0)
    {
        vertex_level = get_vertex_level(parent_id) + 1;
        children[parent_id].push_back(vertex_id);
    }
    else
    {
        vertex_level = 0;
    }

    level.push_back(vertex_level);

    if (num_levels() < vertex_level+1)
        levelset.resize(vertex_level+1);

    levelset[vertex_level].push_back(vertex_id);

    return vertex_id;
}

void CoverTree::set_max_radius()
{
    max_radius = -1;

    for (int64_t i = 0; i < num_points(); ++i)
    {
        max_radius = std::max(max_radius, distance(points[0], points[i]));
    }
}

std::vector<std::tuple<int64_t, std::vector<int64_t>>> CoverTree::init_build_stack()
{
    std::vector<std::tuple<int64_t, std::vector<int64_t>>> stack;
    stack.emplace_back(add_vertex(0, -1), std::vector<int64_t>(num_points()));

    for (int64_t i = 0; i < num_points(); ++i)
        std::get<1>(stack.back())[i] = i;

    return std::move(stack);
}

std::vector<std::tuple<int64_t, std::vector<int64_t>>>
CoverTree::get_next_parents(int64_t parent_id, const std::vector<int64_t>& descendants)
{
    size_t n_descendants = descendants.size();
    std::vector<double> dists(n_descendants);
    std::vector<int64_t> closest(n_descendants);

    for (size_t i = 0; i < n_descendants; ++i)
    {
        dists[i] = point_dist(descendants[0], descendants[i]);
        closest[i] = descendants[0];
    }

    std::vector<int64_t> child_descendants = {get_point_id(parent_id)};
    double parent_ball_radius = vertex_ball_radius(parent_id);

    for (size_t k = 0; k < n_descendants; ++k)
    {
        size_t farthest = std::distance(dists.begin(), std::max_element(dists.begin(), dists.end()));

        if (dists[farthest] <= 0.5 * parent_ball_radius)
            break;

        int64_t next_child_id = descendants[farthest];
        child_descendants.push_back(next_child_id);

        for (size_t j = 0; j < n_descendants; ++j)
        {
            double lastdist = dists[j];
            double curdist = point_dist(descendants[j], next_child_id);

            if (curdist <= lastdist)
            {
                dists[j] = curdist;
                closest[j] = next_child_id;
            }
        }
    }

    std::vector<std::tuple<int64_t, std::vector<int64_t>>> next_parents;

    for (auto itr = child_descendants.begin(); itr != child_descendants.end(); ++itr)
    {
        int64_t child_pt = *itr;
        std::vector<int64_t> next_descendants = {child_pt};

        for (size_t j = 0; j < n_descendants; ++j)
            if (closest[j] == child_pt && descendants[j] != child_pt)
                next_descendants.push_back(descendants[j]);

        int64_t next_point_id;

        next_point_id = add_vertex(child_pt, parent_id);
        next_parents.emplace_back(next_point_id, next_descendants);
    }

    return next_parents;
}

void CoverTree::build_tree()
{
    set_max_radius();
    auto stack = init_build_stack();

    std::vector<double> dists(num_points());
    std::vector<int64_t> closest(num_points());

    while (stack.size() != 0)
    {
        int64_t parent_id = std::get<0>(stack.back());
        const auto descendants = std::get<1>(stack.back());
        stack.pop_back();

        size_t n_descendants = descendants.size();
        assert(n_descendants >= 1 && get_point_id(parent_id) == descendants[0]);

        if (n_descendants == 1)
            continue;

        auto next_parents = get_next_parents(parent_id, descendants);

        for (auto it = next_parents.begin(); it != next_parents.end(); ++it)
            stack.push_back(*it);
    }

    assert(is_full());
    assert(is_nested());
    assert(is_covering());
}

bool CoverTree::is_full() const
{
    std::vector<bool> bits(num_points(), false);

    for (size_t i = 0; i < num_vertices(); ++i)
        bits[get_point_id(i)] = true;

    for (size_t i = 0; i < num_points(); ++i)
        if (!bits[i]) return false;

    return true;
}

bool CoverTree::is_nested() const
{
    std::vector<int64_t> stack = {0};
    std::vector<int64_t> mychildren;
    int64_t u, v;
    size_t m;
    bool found;

    while (stack.size() != 0)
    {
        u = stack.back(); stack.pop_back();
        m = get_child_ids(u, mychildren);

        found = false;
        for (size_t i = 0; i < m; ++i)
        {
            v = mychildren[i];
            stack.push_back(v);

            if (get_point_id(u) == get_point_id(v))
            {
                if (found) { return false; }
                found = true;
            }
        }
    }

    return true;
}

bool CoverTree::is_covering() const
{
    std::vector<int64_t> stack = {0};
    std::vector<int64_t> mychildren;
    int64_t u, v;
    size_t m;
    double ball_radius, d;

    while (stack.size() != 0)
    {
        u = stack.back(); stack.pop_back();
        m = get_child_ids(u, mychildren);
        ball_radius = vertex_ball_radius(u);

        for (size_t i = 0; i < m; ++i)
        {
            v = mychildren[i];
            stack.push_back(v);
            d = point_dist(get_point_id(u), get_point_id(v));

            if (d > ball_radius)
                return false;
        }
    }

    return true;
}

double CoverTree::average_vertex_degree(int64_t level) const
{
    const auto& vs = levelset[level];
    size_t sum = 0;

    for (auto vit = vs.begin(); vit != vs.end(); ++vit)
        sum += children[*vit].size();

    return (sum + 0.0) / (vs.size() + 0.0);
}

double CoverTree::get_neighborhood_graph(double radius, std::vector<std::vector<int64_t>>& nng, bool sort) const
{
    int64_t n = num_points();
    nng.resize(num_points());

    double t = -omp_get_wtime();

    #pragma omp parallel for schedule(dynamic, 100)
    for (int64_t u = 0; u < n; ++u)
    {
        radii_query(points[u], radius, nng[u]);
    }

    t += omp_get_wtime();

    if (sort)
        for (auto it = nng.begin(); it != nng.end(); ++it)
            std::sort(it->begin(), it->end());

    return t;
}

void CoverTree::print_info() const
{
    std::cout << "* number of points: " << num_points() << "\n";
    std::cout << "* dimension: " << points.getdim() << "\n";
    std::cout << "* number of vertices: " << num_vertices() << "\n";
    std::cout << "* number of levels: " << num_levels() << "\n\n";

    for (int64_t i = 0; i < num_levels(); ++i)
    {
        printf("level %ld has %lld vertices of average degree %.3f\n", i, levelset[i].size(), average_vertex_degree(i));
        //std::cout << "level " << i << " has " << levelset[i].size() << " vertices of average degree " << average_vertex_degree(i) << "\n";
    }

    std::cout << std::endl;
}
