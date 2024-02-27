#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <omp.h>

void read_points(std::vector<Point>& points, const char *fname);
void read_truth(std::vector<std::vector<int64_t>>& truth, uint64_t n, const char *fname);
double get_neighborhood_graph(CoverTree& ct, double radius, std::vector<std::vector<int64_t>>& nng);

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s <points.bin> <queries.truth> <radius>\n", argv[0]);
        return 0;
    }

    const char *points_fname = argv[1];
    const char *truth_fname = argv[2];
    double radius = atof(argv[3]);

    double t;
    std::vector<Point> points;
    std::vector<std::vector<int64_t>> truth, nng;

    read_points(points, points_fname);
    read_truth(truth, points.size(), truth_fname);

    CoverTree ct(points);
    t = get_neighborhood_graph(ct, radius, nng);

    auto cmp1 = truth.begin();
    auto cmp2 = nng.begin();

    for (uint64_t i = 0; i < ct.num_points(); ++i)
        assert(*cmp1++ == *cmp2++);

    fprintf(stderr, "Successfully built neighborhood graph in %.5f seconds\n", t);

    return 0;
}

void read_points(std::vector<Point>& points, const char *fname)
{
    points.clear();

    FILE *f = fopen(fname, "rb");
    uint64_t n, d;

    fread(&n, sizeof(uint64_t), 1, f);
    fread(&d, sizeof(uint64_t), 1, f);

    std::vector<double> p(d);

    for (uint64_t i = 0; i < n; ++i)
    {
        fread(p.data(), sizeof(double), d, f);
        points.emplace_back(p);
    }

    fclose(f);
}

void read_truth(std::vector<std::vector<int64_t>>& truth, uint64_t n, const char *fname)
{
    truth.clear(); truth.resize(n);
    std::ifstream f(fname);

    uint64_t m, val;

    for (uint64_t i = 0; i < n; ++i) 
    {
        std::string line;
        std::getline(f, line);
        std::istringstream ss(line);
        ss >> m;

        for (uint64_t j = 0; j < m; ++j)
        {
            ss >> val;
            truth[i].push_back(val);
        }
    }

    f.close();
}

double get_neighborhood_graph(CoverTree& ct, double radius, std::vector<std::vector<int64_t>>& nng)
{
    double t;
    nng.resize(ct.num_points());
    const auto& points = ct.get_points();

    t = -omp_get_wtime();
    for (int64_t u = 0; u < ct.num_points(); ++u)
    {
        ct.radii_query(points[u], radius, nng[u]);
        std::sort(nng[u].begin(), nng[u].end());
    }
    t += omp_get_wtime();

    return t;
}
