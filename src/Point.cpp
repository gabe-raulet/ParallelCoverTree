#include "Point.h"
#include <random>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <omp.h>

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

PointStore PointStore::pgenerate_random_points(size_t n, int d, int seed, double min, double max, int nthreads)
{
    nthreads = std::min(nthreads, omp_get_max_threads());
    omp_set_num_threads(nthreads);

    std::random_device rd;
    std::vector<std::default_random_engine> gens;

    for (int i = 1; i <= nthreads; ++i)
    {
        gens.emplace_back(seed < 0? rd() : seed*i);
    }

    size_t tot = d*n;
    double *p = new double[tot];

    #pragma omp parallel shared(gens, nthreads, tot, p)
    {
        assert(nthreads == omp_get_num_threads());
        std::uniform_real_distribution<> dis(min, max);

        int tid = omp_get_thread_num();

        #pragma omp for
        for (size_t i = 0; i < tot; ++i)
            p[i] = dis(gens[tid]);
    }

    std::shared_ptr<double[]> mem;
    mem.reset(p);

    return PointStore(n, d, mem);
}

bool operator==(const Point& pt1, const Point& pt2)
{
    if (pt1.d != pt2.d) return false;

    double *p1 = pt1.data, *p2 = pt2.data;

    for (int i = 0; i < pt1.d; ++i)
        if (*p1++ != *p2++)
            return false;

    return true;
}

bool operator==(const PointStore& lhs, const PointStore& rhs)
{
    if (lhs.getsize() != rhs.getsize() || lhs.getdim() != rhs.getdim())
        return false;

    for (size_t i = 0; i < lhs.getsize(); ++i)
        if (!(lhs[i] == rhs[i]))
            return false;

    return true;
}
