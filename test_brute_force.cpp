#include "Point.h"
#include "BruteForce.h"
#include "read_args.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
    int num_points = 10;
    int seed = -1;
    int d = 2;
    double radius = 0.25;

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

    num_points = find_int_arg(argc, argv, "-n", num_points);
    d = find_int_arg(argc, argv, "-d", d);
    radius = find_double_arg(argc, argv, "-r", radius);
    seed = find_int_arg(argc, argv, "-s", seed);

    std::vector<Point> points = Point::random_points(num_points, d, seed);

    BruteForce bf(points);
    std::vector<int64_t> ids;

    for (size_t i = 0; i < points.size(); ++i)
    {
        bf.radii_query(points[i], radius, ids);
        std::copy(ids.begin(), ids.end(), std::ostream_iterator<int>(std::cout, ", "));
        std::cout << std::endl;
    }

    return 0;
}

