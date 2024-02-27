#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <cassert>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <omp.h>

double get_neighborhood_graph(CoverTree& ct, double radius, std::vector<std::vector<int64_t>>& nng);
void write_nng(const std::vector<std::vector<int64_t>>& nng, double radius, const char *output);
void write_points(const std::vector<Point>& points, const char *output);

int main(int argc, char *argv[])
{
    int64_t n;
    int d = 2;
    double r0 = 0.05;
    double inc = 0.15;
    int filtsize = 5;
    int seed = -1;
    const char *output = "output";

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -d INT    point dimension [%d]\n", d);
        fprintf(stderr, "         -r FLOAT  start radius [%.4f]\n", r0);
        fprintf(stderr, "         -i FLOAT  radius increment [%.4f]\n", inc);
        fprintf(stderr, "         -f INT    filter size [%d]\n", filtsize);
        fprintf(stderr, "         -s INT    seed [default: random]\n");
        fprintf(stderr, "         -o STR    output prefix [default: '%s']\n", output);
        fprintf(stderr, "         -h        help message\n");
        return 0;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    d = read_int_arg(argc, argv, "-d", &d);
    r0 = read_double_arg(argc, argv, "-r", &r0);
    inc = read_double_arg(argc, argv, "-i", &inc);
    filtsize = read_int_arg(argc, argv, "-f", &filtsize);
    seed = read_int_arg(argc, argv, "-s", &seed);
    output = read_string_arg(argc, argv, "-o", (char**)&output);

    std::vector<Point> points = Point::random_points(n, d, seed);
    write_points(points, output);
    CoverTree ct(points);
    double radius = r0;
    double telapsed;

    for (int step = 0; step < filtsize; step++, radius += inc)
    {
        std::vector<std::vector<int64_t>> nng;
        telapsed = get_neighborhood_graph(ct, radius, nng);
        write_nng(nng, radius, output);
    }

    return 0;
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

void write_nng(const std::vector<std::vector<int64_t>>& nng, double radius, const char *output)
{
    std::string fname = std::string(output) + ".r" + std::to_string(radius);    
    std::ofstream f(fname);

    for (size_t i = 0; i < nng.size(); ++i)
    {
        std::copy(nng[i].begin(), nng[i].end(), std::ostream_iterator<int64_t>(f, " "));
        f << "\n";
    }

    f.close();
}

void write_points(const std::vector<Point>& points, const char *output)
{
    assert(points.size() != 0);

    uint64_t n = points.size();
    uint64_t d = points.back().getdim();
    std::string fname = std::string(output) + ".points.bin";

    FILE *f = fopen(fname.c_str(), "wb");

    fwrite(&n, sizeof(uint64_t), 1, f);
    fwrite(&d, sizeof(uint64_t), 1, f);

    for (uint64_t i = 0; i < n; ++i)
    {
        fwrite(points[i].getdata(), sizeof(double), d, f);
    }

    fclose(f);
}
