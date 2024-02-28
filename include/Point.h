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

    static std::vector<Point> random_points(size_t num_points, int d, int seed = -1);
    static void write_points(const std::vector<Point>& points, const char *fname);
    static void read_points(std::vector<Point>& points, const char *fname);

    friend double distance(const Point& pt1, const Point& pt2);
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
    PointStore(size_t n, int d, std::shared_ptr<double> mem);
    PointStore(const PointStore& rhs);

    int getdim() const { return dim; }
    size_t getsize() const { return points.size(); }

    const Point& operator[](size_t idx) const { return points[idx]; }

    void write_points(const char *fname);
    static PointStore read_points(const char *fname);
    static PointStore generate_random_points(size_t n, int d, int seed, double min, double max);

private:
    int dim;
    std::shared_ptr<double> mem;
    std::vector<Point> points;
};

#endif
