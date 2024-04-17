#include "Point.h"
#include <cmath>

Point::Point() { data[0] = data[1] = 0; }
Point::Point(const float *p) { data[0] = p[0], data[1] = p[1]; }
Point::Point(const Point& rhs) : Point(rhs.data) {}
Point::Point(const pair<float, float>& p) { data[0] = p.first, data[1] = p.second; }

Point& Point::operator=(const Point& rhs)
{
    Point tmp(rhs);
    tmp.swap(*this);
    return *this;
}

void Point::swap(Point& rhs)
{
    ::swap(data, rhs.data);
}

double Point::distance(const Point& rhs) const
{
    double dx = data[0] - rhs.data[0];
    double dy = data[1] - rhs.data[1];
    return sqrt(dx*dx + dy*dy);
}

double distance(const Point& p, const Point& q)
{
    return p.distance(q);
}

vector<Point> Point::random_points(int64_t num_points, double var, int seed)
{
    vector<float> point_data(num_points * Point::dim);
    vector<Point> points;

    random_device rd;
    default_random_engine gen(seed < 0? rd() : 17*seed);
    normal_distribution dis{0.0, sqrt(var)};
    generate(point_data.begin(), point_data.end(), [&]() { return dis(gen); });

    points.reserve(num_points);

    for (int64_t i = 0; i < num_points; ++i)
        points.emplace_back(&point_data[i*Point::dim]);

    return points;
}
