#include "CoverTree.h"
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <random>
#include <tuple>
#include <assert.h>
#include <stdio.h>

double distance(const float *p, const float *q, int d)
{
    double sum = 0., di;

    for (int i = 0; i < d; ++i)
    {
        di = p[i] - q[i];
        sum += di*di;
    }

    return std::sqrt(sum);
}

CoverTree::CoverTree(const float *p, index_t n, int d, double base)
    : base(base), pointmem(new float[d*n]), npoints(n), d(d)
{ 
    std::copy(p, p + n*d, pointmem.get());
    build_tree();
}

bool CoverTree::is_full() const
{
    std::vector<bool> bits(num_points(), false);

    for (size_t i = 0; i < num_vertices(); ++i)
        bits[get_point_id(i)] = true;

    return std::all_of(bits.begin(), bits.end(), [](bool item) { return item == true; });
}

bool CoverTree::is_nested() const
{
    std::vector<index_t> stack = {0};
    std::vector<index_t> mychildren;
    index_t u, v;
    index_t m;
    bool found;

    while (stack.size() != 0)
    {
        u = stack.back(); stack.pop_back();
        mychildren = get_child_ids(u);

        found = false;

        for (index_t v : mychildren)
        {
            stack.push_back(v);

            if (get_point_id(u) == get_point_id(v))
            {
                if (found) return false;
                found = true;
            }
        }
    }

    return true;
}

bool CoverTree::is_covering() const
{
    std::vector<index_t> stack = {0};
    std::vector<index_t> mychildren;
    index_t u, v;
    index_t m;
    double ball_radius, d;

    while (stack.size() != 0)
    {
        u = stack.back(); stack.pop_back();
        mychildren = get_child_ids(u);
        ball_radius = vertex_ball_radius(u);

        for (index_t v : mychildren)
        {
            stack.push_back(v);
            d = point_dist(get_point_id(u), get_point_id(v));

            if (d > ball_radius)
                return false;
        }
    }

    return true;
}

std::vector<index_t> CoverTree::radii_query(const float *query, double radius) const
{
    std::unordered_set<index_t> idset;
    std::vector<index_t> stack = {0};

    index_t u, v;
    std::vector<index_t> mychildren;

    while (stack.size() != 0)
    {
        u = stack.back(); stack.pop_back();
        mychildren = get_child_ids(u);

        for (index_t v : mychildren)
        {
            const float *p = &pointmem[get_point_id(v)*d];

            if (distance(query, p, d) <= radius + max_radius*vertex_ball_radius(v))
                stack.push_back(v);

            if (distance(query, p, d) <= radius)
                idset.insert(get_point_id(v));
        }
    }

    return std::vector<index_t>(idset.begin(), idset.end());
}

void CoverTree::build_tree()
{
    set_max_radius();

    std::vector<std::tuple<index_t, std::vector<index_t>>> stack;
    stack.emplace_back(add_vertex(0, -1), std::vector<index_t>(npoints));

    for (index_t i = 0; i < npoints; ++i)
        std::get<1>(stack.back())[i] = i;

    while (stack.size() != 0)
    {
        index_t parent_id = std::get<0>(stack.back());
        const auto& descendants = std::get<1>(stack.back());
        index_t n_descendants = descendants.size();
        assert(n_descendants >= 1 && get_point_id(parent_id) == descendants[0]);

        if (n_descendants == 1)
        {
            stack.pop_back();
            continue;
        }

        std::vector<double> dists(n_descendants);
        std::vector<index_t> closest(n_descendants, descendants[0]);

        auto it = descendants.begin();
        std::generate(dists.begin(), dists.end(), [&]() { return point_dist(descendants[0], *it++); });

        if (std::all_of(dists.begin(), dists.end(), [](double d) { return d == 0; }))
        {
            for (index_t duplicate_child_pt : descendants)
                add_vertex(duplicate_child_pt, parent_id);

            stack.pop_back();
            continue;
        }

        std::vector<index_t> child_descendants = {get_point_id(parent_id)};
        double parent_ball_radius = vertex_ball_radius(parent_id);

        for (index_t k = 0; k < n_descendants; ++k)
        {
            index_t farthest = std::distance(dists.begin(), std::max_element(dists.begin(), dists.end()));

            if (dists[farthest] <= (parent_ball_radius/base))
                break;

            index_t next_child_id = descendants[farthest];
            child_descendants.push_back(next_child_id);

            for (index_t j = 0; j < n_descendants; ++j)
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

        std::vector<std::tuple<index_t, std::vector<index_t>>> next_parents;

        for (index_t child_pt : child_descendants)
        {
            std::vector<index_t> next_descendants = {child_pt};

            for (index_t j = 0; j < n_descendants; ++j)
                if (closest[j] == child_pt && descendants[j] != child_pt)
                    next_descendants.push_back(descendants[j]);

            index_t next_point_id = add_vertex(child_pt, parent_id);
            next_parents.emplace_back(next_point_id, next_descendants);
        }

        stack.pop_back();
        stack.reserve(stack.size() + next_parents.size());
        std::copy(next_parents.begin(), next_parents.end(), std::back_inserter(stack));
    }
}

index_t CoverTree::add_vertex(index_t point_id, index_t parent_id)
{
    index_t vertex_level;
    index_t vertex_id = num_vertices();

    pt.push_back(point_id);
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
    levelset.resize(std::max(vertex_level+1, static_cast<index_t>(levelset.size())));
    levelset[vertex_level].push_back(vertex_id);

    return vertex_id;
}

void CoverTree::set_max_radius()
{
    max_radius = -1;

    for (index_t i = 0; i < npoints; ++i)
    {
        max_radius = std::max(max_radius, distance(&pointmem[0], &pointmem[i*d], d));
    }
}

double CoverTree::point_dist(index_t id1, index_t id2) const
{
    return distance(&pointmem[d*id1], &pointmem[d*id2], d) / max_radius;
}

double CoverTree::vertex_ball_radius(index_t vertex_id) const
{
    return std::pow(base, -1. * get_vertex_level(vertex_id));
}

void CoverTree::print_info() const
{
    std::cout << "* number of points: " << num_points() << "\n";
    std::cout << "* dimension: " << getdim() << "\n";
    std::cout << "* number of vertices: " << num_vertices() << "\n";
    std::cout << "* number of levels: " << num_levels() << "\n\n";

    for (int64_t i = 0; i < num_levels(); ++i)
    {
        double average_vertex_degree = 0;

        for (const auto& v : levelset[i])
        {
            average_vertex_degree += children[v].size();
        }

        average_vertex_degree /= levelset[i].size();
        printf("level %ld has %lld vertices of average degree %.3f\n", i, levelset[i].size(), average_vertex_degree);
    }

    std::cout << std::endl;

}

void CoverTree::read_from_file(const char *fname)
{
    FILE *f = fopen(fname, "r");
    assert(f != NULL);

    index_t n_vertices;

    fread(&npoints, sizeof(index_t), 1, f); assert(npoints >= 1);
    fread(&d, sizeof(int), 1, f); assert(d >= 1);

    pointmem.reset(new float[d*npoints]);
    fread(pointmem.get(), sizeof(float), d*npoints, f);
    fread(&max_radius, sizeof(double), 1, f);
    fread(&base, sizeof(double), 1, f);
    fread(&n_vertices, sizeof(index_t), 1, f);

    std::vector<index_t> parent(n_vertices);
    pt.resize(n_vertices);
    level.resize(n_vertices);

    fread(pt.data(), sizeof(index_t), n_vertices, f);
    fread(parent.data(), sizeof(index_t), n_vertices, f);
    fread(level.data(), sizeof(index_t), n_vertices, f);

    index_t num_levels = std::accumulate(level.begin(), level.end(), static_cast<index_t>(0), [](auto a, auto b) { return std::max(a,b); }) + 1;

    children.resize(n_vertices);
    levelset.resize(num_levels);
    levelset[level[0]].push_back(0);

    for (index_t v = 1; v < n_vertices; ++v)
    {
        children[parent[v]].push_back(v);
        levelset[level[v]].push_back(v);
    }

    fclose(f);
}

void CoverTree::write_to_file(const char *fname) const
{
    index_t n_vertices = num_vertices();
    std::vector<index_t> parent(n_vertices);
    parent[0] = -1;

    for (index_t i = 0; i < n_vertices; ++i)
    {
        const auto& child_ids = children[i];

        for (index_t child_id : child_ids)
            parent[child_id] = i;
    }

    FILE *f = fopen(fname, "w");

    fwrite(&npoints, sizeof(index_t), 1, f);
    fwrite(&d, sizeof(int), 1, f);
    fwrite(pointmem.get(), sizeof(float), d*npoints, f);
    fwrite(&max_radius, sizeof(double), 1, f);
    fwrite(&base, sizeof(double), 1, f);
    fwrite(&n_vertices, sizeof(index_t), 1, f);
    fwrite(pt.data(), sizeof(index_t), n_vertices, f);
    fwrite(parent.data(), sizeof(index_t), n_vertices, f);
    fwrite(level.data(), sizeof(index_t), n_vertices, f);

    fclose(f);
}
