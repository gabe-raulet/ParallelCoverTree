#include "Point.h"
#include <random>
#include <string>
#include <assert.h>
#include <stdio.h>

Point::Point(int d) : isowner(true), d(d), data(new double[d]()) {}

Point::Point(const double *p, int d, bool isowner) : d(d), isowner(isowner)
{
    if (isowner)
    {
        data = new double[d];
        std::copy(p, p + d, data);
    }
    else
    {
        data = const_cast<double*>(p);
    }
}

Point::Point(const Point& rhs) : Point(rhs.data, rhs.d) {}

//std::vector<Point> Point::random_points(size_t num_points, int d, int seed)
//{
    //std::random_device rd;
    //std::default_random_engine gen(seed < 0? rd() : seed);
    //std::uniform_real_distribution<> dis(-1., 1.);

    //double *mem = new double[d*num_points];

    //for (size_t i = 0; i < d*num_points; ++i)
        //mem[i] = dis(gen);

    //std::vector<Point> points;
    //const double *p = mem;

    //for (size_t i = 0; i < num_points; ++i)
    //{
        //points.emplace_back(p, d);
        //p += d;
    //}

    //delete[] mem;

    //return points;
//}

//void Point::write_points(const std::vector<Point>& points, const char *fname)
//{
    //assert(points.size() != 0);

    //uint64_t n = points.size();
    //uint64_t d = points.back().getdim();

    //FILE *f = fopen(fname, "wb");

    //fwrite(&n, sizeof(uint64_t), 1, f);
    //fwrite(&d, sizeof(uint64_t), 1, f);

    //for (uint64_t i = 0; i < n; ++i)
    //{
        //fwrite(points[i].getdata(), sizeof(double), d, f);
    //}

    //fclose(f);
//}

//void Point::read_points(std::vector<Point>& points, const char *fname)
//{
    //points.clear();

    //FILE *f = fopen(fname, "rb");
    //uint64_t n, d;

    //fread(&n, sizeof(uint64_t), 1, f);
    //fread(&d, sizeof(uint64_t), 1, f);

    //std::vector<double> p(d);

    //for (uint64_t i = 0; i < n; ++i)
    //{
        //fread(p.data(), sizeof(double), d, f);
        //points.emplace_back(p.data(), d);
    //}

    //fclose(f);
//}

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

PointStore::PointStore(size_t n, int d, std::shared_ptr<double[]> mem) : dim(d), mem(mem)
{
    points.reserve(n);
    double *p = mem.get();

    for (size_t i = 0; i < n; ++i)
    {
        points.emplace_back(p, dim, false);
        p += dim;
    }
}

PointStore::PointStore(const PointStore& rhs) : dim(rhs.dim), mem(rhs.mem), points(rhs.points) {}

void PointStore::write_points(const char *fname)
{
    FILE *f = fopen(fname, "wb");

    uint64_t n = getsize();
    uint64_t d = getdim();

    fwrite(&n, sizeof(uint64_t), 1, f);
    fwrite(&d, sizeof(uint64_t), 1, f);
    fwrite(mem.get(), sizeof(double), n*d, f);

    fclose(f);
}

PointStore PointStore::read_points(const char *fname)
{
    FILE *f = fopen(fname, "rb");
    assert(f != NULL);

    uint64_t n, d;

    fread(&n, sizeof(uint64_t), 1, f);
    fread(&d, sizeof(uint64_t), 1, f);

    std::shared_ptr<double[]> mem;
    mem.reset(new double[n*d]);

    fread(mem.get(), sizeof(double), n*d, f);
    fclose(f);

    return PointStore(n, d, mem);
}

PointStore PointStore::generate_random_points(size_t n, int d, int seed, double min, double max)
{
    std::random_device rd;
    std::default_random_engine gen(seed < 0? rd() : seed);
    std::uniform_real_distribution<> dis(min, max);

    std::shared_ptr<double[]> mem;
    mem.reset(new double[n*d]);
    double *p = mem.get();

    for (size_t i = 0; i < d*n; ++i)
        p[i] = dis(gen);

    return PointStore(n, d, mem);
}
