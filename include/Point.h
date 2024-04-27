#ifndef POINT_H_
#define POINT_H_

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <mpi.h>

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

    void fill_dest(float *dest) const { dest[0] = data[0], dest[1] = data[1]; }

    string repr() const;

    static void create_mpi_dtype(MPI_Datatype *MPI_POINT);

    static vector<Point> random_points(int64_t num_points, double var, int seed);
    static vector<Point> dist_random_points(int64_t num_points, double var, int seed, int root, MPI_Comm comm);

    static vector<Point> scatter(const vector<Point>& points, int root, MPI_Comm comm);

    static vector<Point> from_file(const char *fname);
    static void to_file(const vector<Point>& points, const char *fname);

    friend ostream& operator<<(ostream& os, const Point& pt)
    {
        os << setprecision(3) << "[" << pt.data[0] << ", " << pt.data[1] << "]";
        return os;
    }

private:
    float data[dim];
};

#endif
