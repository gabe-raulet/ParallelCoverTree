#include "Point.h"
#include <cmath>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <numeric>

Point::Point() : data{} {}
Point::Point(const float *p) { copy(p, p + dim, data); }
Point::Point(const Point& rhs) : Point(rhs.data) {}

void Point::fill_dest(float *dest) const { copy(data, data + dim, dest); }

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
    double sum = 0.0, delta;

    for (int i = 0; i < dim; ++i)
    {
        delta = data[i] - rhs.data[i];
        sum += delta *delta;
    }

    return sqrt(sum);
}

double distance(const Point& p, const Point& q)
{
    return p.distance(q);
}

vector<Point> Point::random_points(int64_t num_points, double var, int seed)
{
    vector<float> point_data(num_points * dim);
    vector<Point> points;

    default_random_engine gen(17*seed);
    normal_distribution dis{0.0, sqrt(var)};
    generate(point_data.begin(), point_data.end(), [&]() { return dis(gen); });

    points.reserve(num_points);

    for (int64_t i = 0; i < num_points; ++i)
        points.emplace_back(&point_data[i*dim]);

    return points;
}

void Point::create_mpi_dtype(MPI_Datatype *MPI_POINT)
{
    MPI_Type_contiguous(dim, MPI_FLOAT, MPI_POINT);
    MPI_Type_commit(MPI_POINT);
}

vector<Point> Point::scatter(const vector<Point>& points, int root, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    vector<Point> mypoints;
    vector<int> sendcounts, displs;
    int recvcount;

    if (myrank == root)
    {
        sendcounts.resize(nprocs);
        displs.resize(nprocs);
        fill(sendcounts.begin(), sendcounts.end(), points.size()/nprocs);
        sendcounts.back() = points.size() - (nprocs-1)*sendcounts.back();
        displs.front() = 0;
        partial_sum(sendcounts.begin(), sendcounts.end()-1, displs.begin()+1);
    }

    MPI_Scatter(sendcounts.data(), 1, MPI_INT, &recvcount, 1, MPI_INT, root, comm);

    MPI_Datatype MPI_POINT;
    create_mpi_dtype(&MPI_POINT);

    mypoints.resize(recvcount);
    MPI_Scatterv(points.data(), sendcounts.data(), displs.data(), MPI_POINT, mypoints.data(), recvcount, MPI_POINT, root, comm);

    MPI_Type_free(&MPI_POINT);
    return mypoints;
}

vector<Point> Point::dist_random_points(int64_t num_points, double var, int seed, int root, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    vector<Point> points;

    if (myrank == root) points = random_points(num_points, var, seed);

    return scatter(points, root, comm);
}

vector<Point> Point::from_file(const char *fname)
{
    int64_t n;
    FILE *f;
    vector<float> point_data;
    vector<Point> points;

    f = fopen(fname, "r");
    assert(f != NULL);

    fread(&n, sizeof(int64_t), 1, f);
    point_data.resize(n*dim);
    fread(point_data.data(), sizeof(float), dim*n, f);
    fclose(f);

    points.reserve(n);

    for (int64_t i = 0; i < n; ++i)
    {
        points.emplace_back(&point_data[i*dim]);
    }

    return points;
}

void Point::to_file(const vector<Point>& points, const char *fname)
{
    int64_t n;
    FILE *f;
    vector<float> point_data;
    float *dest;

    n = points.size();
    point_data.resize(n*dim);
    dest = point_data.data();

    for (int64_t i = 0; i < n; ++i)
    {
        const Point& p = points[i];
        p.fill_dest(dest);
        dest += dim;
    }

    f = fopen(fname, "w");
    assert(f != NULL);

    fwrite(&n, sizeof(int64_t), 1, f);
    fwrite(point_data.data(), sizeof(float), dim*n, f);

    fclose(f);
}
