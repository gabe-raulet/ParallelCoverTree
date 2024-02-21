#ifndef POINT_H_
#define POINT_H_

#include <vector>
#include <iostream>

class Point
{
public:
    Point(int d);
    Point(const double *p, int d);
    Point(const std::vector<double>& p);
    Point(const Point& rhs);

    static std::vector<Point> random_points(size_t num_points, int d, int seed = -1);

    friend double distance(const Point& pt1, const Point& pt2);
    friend std::ostream& operator<<(std::ostream& stream, const Point& p);

private:
    double *data;
    int d;
};

#endif
