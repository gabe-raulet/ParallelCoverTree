#include "CoverTree.h"
#include <list>
#include <algorithm>
#include <unordered_set>
#include <iterator>
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

std::vector<std::vector<index_t>> greedy_partitioning(const std::vector<index_t>& tasks, int nprocs)
{
    size_t num_tasks = tasks.size();
    std::vector<index_t> loads(nprocs, 0);
    std::vector<std::vector<index_t>> partitions(nprocs);

    for (index_t i = 0; i < num_tasks; ++i)
    {
        int proc = std::distance(loads.begin(), std::min_element(loads.begin(), loads.end()));
        loads[proc] += tasks[proc];
        partitions[proc].push_back(i);
    }

    return partitions;
}

CoverTree::CoverTree(const float *p, index_t n, int d, double base)
    : base(base), pointmem(new float[d*n]), npoints(n), d(d)
{
    std::copy(p, p + n*d, pointmem.get());
    build_tree();
}

bool CoverTree::is_full() const
{
    std::vector<bool> bits(npoints, false);

    for (size_t i = 0; i < num_vertices(); ++i)
        bits[pt[i]] = true;

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
        mychildren = children[u];

        found = false;

        for (index_t v : mychildren)
        {
            stack.push_back(v);

            if (pt[u] == pt[v])
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
        mychildren = children[u];
        ball_radius = vertex_ball_radius(u);

        for (index_t v : mychildren)
        {
            stack.push_back(v);
            d = point_dist(pt[u], pt[v]);

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
        mychildren = children[u];

        for (index_t v : mychildren)
        {
            const float *p = &pointmem[pt[v]*d];

            if (distance(query, p, d) <= radius + max_radius*vertex_ball_radius(v))
                stack.push_back(v);

            if (distance(query, p, d) <= radius)
                idset.insert(pt[v]);
        }
    }

    return std::vector<index_t>(idset.begin(), idset.end());
}

std::vector<std::tuple<index_t, std::vector<index_t>>>
CoverTree::compute_child_points(index_t parent_id, const std::vector<index_t>& descendants) const
{
    std::vector<std::tuple<index_t, std::vector<index_t>>> child_info;
    index_t n_descendants = descendants.size();
    assert(n_descendants >= 1 && pt[parent_id] == descendants[0]);

    if (n_descendants == 1) return child_info;

    std::vector<double> dists(n_descendants);
    std::vector<index_t> closest(n_descendants, descendants[0]);

    auto it = descendants.begin();
    std::generate(dists.begin(), dists.end(), [&]() { return point_dist(descendants[0], *it++); });

    if (std::all_of(dists.begin(), dists.end(), [](double d) { return d == 0; }))
    {
        for (index_t duplicate_child_pt : descendants)
        {
            child_info.emplace_back(duplicate_child_pt, std::vector<index_t>());
        }
        return child_info;
    }

    std::vector<index_t> child_descendants = {pt[parent_id]};
    double parent_ball_radius = vertex_ball_radius(parent_id);

    for (index_t k = 0; k < n_descendants; ++k)
    {
        index_t farthest = std::distance(dists.begin(), std::max_element(dists.begin(), dists.end()));

        if (dists[farthest] <= (parent_ball_radius / base))
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

    for (index_t child_pt : child_descendants)
    {
        std::vector<index_t> next_descendants = {child_pt};

        for (index_t j = 0; j < n_descendants; ++j)
            if (closest[j] == child_pt && descendants[j] != child_pt)
                next_descendants.push_back(descendants[j]);

        child_info.emplace_back(child_pt, next_descendants);
    }

    return child_info;
}

void CoverTree::build_tree()
{
    set_max_radius();

    std::list<std::tuple<index_t, std::vector<index_t>>> queue;

    index_t root_id = add_vertex(0, -1);
    std::vector<index_t> root_descendants(npoints);
    std::iota(root_descendants.begin(), root_descendants.end(), 0);
    queue.emplace_back(root_id, root_descendants);

    std::vector<std::vector<index_t>> descendant_counts;

    while (queue.size() != 0)
    {
        index_t parent_id = std::get<0>(queue.front());
        const auto& descendants = std::get<1>(queue.front());
        const auto& child_info = compute_child_points(parent_id, descendants);

        descendant_counts.resize(level[parent_id]+1);
        descendant_counts[level[parent_id]].push_back(descendants.size());

        queue.pop_front();

        for (const auto& item : child_info)
        {
            const auto& child_pt = std::get<0>(item);
            const auto& next_descendants = std::get<1>(item);

            if (next_descendants.size() == 0)
                continue;

            index_t next_vertex_id = add_vertex(child_pt, parent_id);
            queue.emplace_back(next_vertex_id, next_descendants);
        }
    }

    print_info(stderr);
    std::cerr << std::endl;

    auto level_set = get_level_set();
    std::vector<index_t> num_leaves;

    for (index_t l = 0; l < descendant_counts.size(); ++l)
    {
        index_t c = 0;

        for (auto v : level_set[l])
        {
            if (children[v].size() == 0)
                c++;
        }

        num_leaves.push_back(c);
    }

    for (index_t l = 0; l < descendant_counts.size(); ++l)
    {
        index_t total = std::accumulate(descendant_counts[l].begin(), descendant_counts[l].end(), static_cast<index_t>(0), [](auto a, auto b) { return a + b; });
        std::cerr << "level=" << l << " :: total=" << total << " :: leaves=" << num_leaves[l] << " :: ";
        std::copy(descendant_counts[l].begin(), descendant_counts[l].end(), std::ostream_iterator<index_t>(std::cerr, " "));
        std::cerr << std::endl;
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
    return std::pow(base, -1. * level[vertex_id]);
}

std::vector<std::vector<index_t>> CoverTree::get_level_set() const
{
    index_t num_levels = std::accumulate(level.begin(), level.end(), static_cast<index_t>(0), [](auto a, auto b) { return std::max(a,b); }) + 1;
    std::vector<std::vector<index_t>> level_set(num_levels);

    for (index_t i = 0; i < pt.size(); ++i)
        level_set[level[i]].push_back(i);

    return level_set;
}

void CoverTree::print_info(FILE *f) const
{
    const auto& level_set = get_level_set();
    index_t num_levels = level_set.size();

    fprintf(f, "* number of points: %lld\n", npoints);
    fprintf(f, "* dimension: %d\n", d);
    fprintf(f, "* number of vertices: %lu\n", pt.size());
    fprintf(f, "* number of levels: %lld\n\n", num_levels);

    for (index_t i = 0; i < num_levels; ++i)
    {
        double average_vertex_degree = 0;

        for (const auto& v : level_set[i])
        {
            average_vertex_degree += children[v].size();
        }

        average_vertex_degree /= level_set[i].size();
        fprintf(f, "level %lld has %lu vertices of average degree %.3f\n", i, level_set[i].size(), average_vertex_degree);
    }

    fprintf(f, "\n");
    fflush(f);
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

    for (index_t v = 1; v < n_vertices; ++v)
        children[parent[v]].push_back(v);

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
