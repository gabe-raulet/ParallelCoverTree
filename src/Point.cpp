#include "Point.h"
#include <random>
#include <string>
#include <assert.h>
#include <stdio.h>

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

void Point::write_points(const std::vector<Point>& points, const char *fname)
{
    assert(points.size() != 0);

    uint64_t n = points.size();
    uint64_t d = points.back().getdim();

    FILE *f = fopen(fname, "wb");

    fwrite(&n, sizeof(uint64_t), 1, f);
    fwrite(&d, sizeof(uint64_t), 1, f);

    for (uint64_t i = 0; i < n; ++i)
    {
        fwrite(points[i].getdata(), sizeof(double), d, f);
    }

    fclose(f);
}

void Point::read_points(std::vector<Point>& points, const char *fname)
{
    points.clear();

    FILE *f = fopen(fname, "rb");
    uint64_t n, d;

    fread(&n, sizeof(uint64_t), 1, f);
    fread(&d, sizeof(uint64_t), 1, f);

    std::vector<double> p(d);

    for (uint64_t i = 0; i < n; ++i)
    {
        fread(p.data(), sizeof(double), d, f);
        points.emplace_back(p);
    }

    fclose(f);
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
