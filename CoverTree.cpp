#include "CoverTree.h"
#include <cmath>

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
    return 0;
}
