#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <algorithm>
#include <cassert>
#include <iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>

int write_points(const std::vector<Point>& points, const char *fname);
int write_queries(CoverTree& ct, double radius, const char *fname);

int main(int argc, char *argv[])
{
    int num_points = 10;
    int seed = -1;
    int d = 2;
    double radius = 0.25;
    char *outprefix = NULL;

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
        std::cerr << "Options: -n INT    number of points [default: " << num_points << "]" << std::endl;
        std::cerr << "         -d INT    point dimension [default: " << d << "]" << std::endl;
        std::cerr << "         -r FLOAT  query radius [default: " << radius << "]" << std::endl;
        std::cerr << "         -s INT    seed [default: random]" << std::endl;
        std::cerr << "         -o STR    output prefix [required]" << std::endl;
        std::cerr << "         -h        help message" << std::endl;
        return -1;
    }

    num_points = find_int_arg(argc, argv, "-n", num_points);
    d = find_int_arg(argc, argv, "-d", d);
    radius = find_double_arg(argc, argv, "-r", radius);
    seed = find_int_arg(argc, argv, "-s", seed);
    outprefix = find_string_arg(argc, argv, "-o", outprefix);

    if (!outprefix)
    {
        std::cerr << "error: missing -o argument!" << std::endl;
        return -1;
    }

    CoverTree ct(Point::random_points(num_points, d, seed));
    const std::vector<Point>& points = ct.get_points();
    std::vector<int64_t> ids;

    std::string pts_fname = std::string(outprefix) + ".points.bin";
    std::string qs_fname = std::string(outprefix) + ".queries.txt";

    write_points(points, pts_fname.c_str());
    write_queries(ct, radius, qs_fname.c_str());

    return 0;
}

int write_points(const std::vector<Point>& points, const char *fname)
{
    assert(points.size() != 0);

    FILE *f = fopen(fname, "wb");

    uint64_t num_points = points.size();
    uint64_t dim = points.back().getdim();

    fwrite(&num_points, sizeof(num_points), 1, f);
    fwrite(&dim, sizeof(dim), 1, f);

    for (size_t i = 0; i < num_points; ++i)
    {
        fwrite(points[i].getdata(), sizeof(double), dim, f);
    }

    fclose(f);
    return 0;
}

int write_queries(CoverTree& ct, double radius, const char *fname)
{
    const std::vector<Point>& points = ct.get_points();

    std::ofstream f(fname);
    std::vector<int64_t> ids;

    for (auto itr = points.begin(); itr != points.end(); ++itr)
    {
        ct.radii_query(*itr, radius, ids);
        std::sort(ids.begin(), ids.end());
        std::copy(ids.begin(), ids.end(), std::ostream_iterator<int>(f, ","));
        f << std::endl;
    }

    f.close();
    return 0;
}
