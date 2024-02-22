#include "BruteForce.h"

size_t BruteForce::radii_query(const Point& query, double radius, std::vector<int64_t>& ids)
{
    ids.clear();

    for (int64_t i = 0; i < num_points(); ++i)
    {
        if (distance(query, points[i]) <= radius)
            ids.push_back(i);
    }

    return ids.size();
}
