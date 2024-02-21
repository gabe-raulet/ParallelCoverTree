#ifndef BRUTE_FORCE_H_
#define BRUTE_FORCE_H_

#include "Point.h"
#include <vector>

class BruteForce
{
public:
    BruteForce(const std::vector<Point>& points) : points(points) {}
    ~BruteForce() {}

    size_t radii_query(const Point& query, double radius, std::vector<int64_t>& ids);
    size_t num_points() const { return points.size(); }

private:
    std::vector<Point> points;
};

#endif
