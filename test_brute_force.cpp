#include "Point.h"
#include "BruteForce.h"
#include "read_args.h"
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
    int num_points = 10;
    int seed = -1;

    if (find_arg_idx(argc, argv, "-h") >= 0)
    {
        std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
        std::cerr << "Options: -n INT   number of points [default: " << num_points << "]" << std::endl;
        std::cerr << "         -s INT   seed [default: random]" << std::endl;
        std::cerr << "         -h       help message" << std::endl;
        return -1;
    }

    num_points = find_int_arg(argc, argv, "-n", num_points);
    seed = find_int_arg(argc, argv, "-s", seed);

    std::cout << "num_points=" << num_points << std::endl;
    std::cout << "seed=" << seed << std::endl;

    return 0;
}

