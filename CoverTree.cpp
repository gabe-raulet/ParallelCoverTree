#include "CoverTree.h"
#include <cmath>
#include <unordered_set>

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
    radius /= max_radius;
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

            if (distance(query, points[get_point_id(v)]) <= radius + vertex_ball_radius(v))
            {
                stack.push_back(v);
            }

            if (distance(query, points[get_point_id(v)]) <= radius)
            {
                idset.insert(get_point_id(v));
            }
        }
    }

    ids.clear();
    ids.assign(idset.begin(), idset.end());
    std::sort(ids.begin(), ids.end());

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
