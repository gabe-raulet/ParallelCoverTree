#ifndef POINT_H_
#define POINT_H_

#include <vector>
#include <iostream>
#include <memory>

class Point
{
public:
    Point(int d);
    Point(const double *p, int d, bool isowner = true);
    Point(const Point& rhs);

    int getdim() const { return d; }
    const void* getdata() const { return static_cast<const void*>(data); }

    friend double distance(const Point& pt1, const Point& pt2);
    friend bool operator==(const Point& pt1, const Point& pt2);
    friend std::ostream& operator<<(std::ostream& stream, const Point& p);

    ~Point() { if (isowner) delete[] data; }

private:
    int d;
    bool isowner;
    double *data;
};

class PointStore
{
public:
    PointStore(size_t n, int d, std::shared_ptr<double[]> mem);
    PointStore(const PointStore& rhs);

    int getdim() const { return dim; }
    size_t getsize() const { return points.size(); }
    const std::vector<Point>& getvector() const { return points; }

    const Point& operator[](size_t idx) const { return points[idx]; }

    void write_points(const char *fname);
    static PointStore read_points(const char *fname);
    static PointStore generate_random_points(size_t n, int d, int seed, double min, double max);
    static PointStore pgenerate_random_points(size_t n, int d, int seed, double min, double max, int nthreads);

    friend bool operator==(const PointStore& lhs, const PointStore& rhs);

private:
    int dim;
    std::shared_ptr<double[]> mem;
    std::vector<Point> points;
};

#endif
