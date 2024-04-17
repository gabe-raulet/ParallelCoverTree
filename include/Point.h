#ifndef POINT_H_
#define POINT_H_

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>

using namespace std;

class Point
{
public:
    static constexpr int dim = 2;

    Point();
    Point(const float *p);
    Point(const Point& rhs);
    Point(const pair<float, float>& p);

    Point& operator=(const Point& rhs);

    void swap(Point& rhs);

    double distance(const Point& rhs) const;
    friend double distance(const Point& p, const Point& q);

    string repr() const;

    static vector<Point> random_points(int64_t num_points, double var, int seed);

    friend ostream& operator<<(ostream& os, const Point& pt)
    {
        os << setprecision(3) << "[" << pt.data[0] << ", " << pt.data[1] << "]";
        return os;
    }

private:
    float data[dim];
};

#endif
