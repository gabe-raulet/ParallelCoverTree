#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <iterator>
#include <iostream>
#include <vector>
#include <omp.h>

int main(int argc, char *argv[])
{
    int64_t num_points = 1024*64;
    int seed = -1;
    int d = 16;
    double radius = 0.5;

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
        std::cerr << "Options: -n INT    number of points [default: " << num_points << "]" << std::endl;
        std::cerr << "         -d INT    point dimension [default: " << d << "]" << std::endl;
        std::cerr << "         -r FLOAT  query radius [default: " << radius << "]" << std::endl;
        std::cerr << "         -s INT    seed [default: random]" << std::endl;
        std::cerr << "         -h        help message" << std::endl;
        return -1;
    }

    num_points = read_formatted_int_arg(argc, argv, "-n", &num_points);
    d = read_int_arg(argc, argv, "-d", &d);
    radius = read_double_arg(argc, argv, "-r", &radius);
    seed = read_int_arg(argc, argv, "-s", &seed);


    std::vector<Point> points = Point::random_points(num_points, d, seed);

    double t;

    t = -omp_get_wtime();
    CoverTree ct(points);
    t += omp_get_wtime();

    std::cout << "Tree construction took " << t << " seconds on " << num_points << " points of dimension " << d << std::endl;

    std::cout << "num_vertices: " << ct.num_vertices() << std::endl;
    std::cout << "num_levels: " << ct.num_levels() << std::endl;
    std::cout << "vertices_per_level: " << ct.vertices_per_level() << std::endl;

    std::vector<int64_t> dummy;

    for (int64_t i = 0; i < ct.num_levels(); ++i)
    {
        std::cout << "level " << i << " has " << ct.get_level_ids(i, dummy) << " vertices" << std::endl;
    }

    return 0;
}

