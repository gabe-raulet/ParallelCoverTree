#include "Point.h"
#include <random>
#include <cassert>

Point::Point(int d) : d(d), data(new double[d]()) {}

Point::Point(const double *p, int d) : d(d), data(new double[d])
{
    std::copy(p, p + d, data);
}

Point::Point(const std::vector<double>& p) : d(p.size()), data(new double[p.size()])
{
    std::copy(p.begin(), p.end(), data);
}

Point::Point(const Point& rhs) : Point(rhs.data, rhs.d) {}

std::vector<Point> Point::random_points(size_t num_points, int d, int seed)
{
    std::random_device rd;
    std::default_random_engine gen(seed < 0? rd() : seed);
    std::uniform_real_distribution<> dis(-1., 1.);

    double *mem = new double[d*num_points];

    for (size_t i = 0; i < d*num_points; ++i)
        mem[i] = dis(gen);

    std::vector<Point> points;
    const double *p = mem;

    for (size_t i = 0; i < num_points; ++i)
    {
        points.emplace_back(p, d);
        p += d;
    }

    delete[] mem;

    return points;
}

double distance(const Point& pt1, const Point& pt2)
{
    assert(pt1.d == pt2.d);

    double sum = 0., di;

    for (int i = 0; i < pt1.d; ++i)
    {
        di = pt1.data[i] - pt2.data[i];
        sum += di*di;
    }

    return std::sqrt(sum);
}

std::ostream& operator<<(std::ostream& stream, const Point& p)
{
    stream << "Point(";
    for (int i = 0; i < p.d-1; ++i)
        stream << p.data[i] << ",";
    stream << p.data[p.d-1] << ")";
    return stream;
}
