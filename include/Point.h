#ifndef POINT_H_
#define POINT_H_

#include <vector>
#include <iostream>

class Point
{
public:
    Point(int d);
    Point(const double *p, int d);
    Point(const Point& rhs);

    int getdim() const { return d; }
    const void* getdata() const { return static_cast<const void*>(data); }

    static std::vector<Point> random_points(size_t num_points, int d, int seed = -1);
    static void write_points(const std::vector<Point>& points, const char *fname);
    static void read_points(std::vector<Point>& points, const char *fname);

    friend double distance(const Point& pt1, const Point& pt2);
    friend std::ostream& operator<<(std::ostream& stream, const Point& p);

private:
    double *data;
    int d;
};

#endif
