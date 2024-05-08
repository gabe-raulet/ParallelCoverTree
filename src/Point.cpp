#include "Point.h"
#include <cmath>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <filesystem>

namespace fs = filesystem;

Point::Point() : data{} {}
Point::Point(const float *p) { copy(p, p + dim, data); }
Point::Point(const Point& rhs) : Point(rhs.data) {}

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
    size_t filesize, n;
    int d;
    FILE *f;
    fs::path path;
    vector<Point> points;
    float buf[dim];

    path = fname;
    filesize = fs::file_size(path);

    f = fopen(fname, "rb");
    fread(&d, sizeof(int), 1, f);
    assert(d == dim);
    fseek(f, SEEK_SET, 0);
    n = filesize / (4*(dim+1));

    points.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        fread(&d, sizeof(int), 1, f);
        fread(buf, sizeof(float), dim, f);
        points.emplace_back(buf);
    }

    fclose(f);
    return points;
}

void Point::to_file(const vector<Point>& points, const char *fname)
{
    size_t n;
    FILE *f;

    n = points.size();
    f = fopen(fname, "wb");

    for (size_t i = 0; i < n; ++i)
    {
        fwrite(&dim, sizeof(int), 1, f);
        fwrite(points[i].getdata(), sizeof(float), dim, f);
    }

    fclose(f);
}
