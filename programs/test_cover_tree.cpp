#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <iterator>
#include <iostream>
#include <vector>
#include <omp.h>

int main(int argc, char *argv[])
{
    const static std::vector<int> dims = {2,4,8};
    int seed = -1; seed = read_int_arg(argc, argv, "-s", &seed);

    for (auto ditr = dims.begin(); ditr != dims.end(); ++ditr)
    {
        int d = *ditr;

        for (int exp = 12; exp <= 22; ++exp)
        {
            int64_t num_points = 1LL << exp;
            std::vector<Point> points = Point::random_points(num_points, d, seed*exp);

            double t;

            t = -omp_get_wtime();
            CoverTree ct = CoverTree::build(points);
            t += omp_get_wtime();

            std::cout << "** number of points: " << num_points << "\n";
            std::cout << "** dimension: " << d << "\n";
            std::cout << "** construction time (seconds): " << t << "\n";
            std::cout << "** inserts per millsecond: " << static_cast<double>(num_points) / (t*1000) << "\n";
            std::cout << "** number of vertices: " << ct.num_vertices() << "\n";
            std::cout << "** number of levels: " << ct.num_levels() << "\n";

            std::vector<int64_t> dummy;

            for (int64_t i = 0; i < ct.num_levels(); ++i)
            {
                std::cout << "level " << i << " has " << ct.get_level_ids(i, dummy) << " vertices of average degree " << ct.average_vertex_degree(i) << "\n";
            }

            std::cout << "\n\n";
            std::cout << std::flush;
        }
    }

    return 0;
}
