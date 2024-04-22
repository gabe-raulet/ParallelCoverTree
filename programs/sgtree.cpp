#include "Point.h"
#include "SGTree.h"
#include "MyTimer.h"
#include "read_args.h"
#include "version.h"
#include <iostream>
#include <chrono>
#include <assert.h>
#include <stdio.h>

using namespace std;

string return_current_time_and_date();
string program_str(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    int64_t n;
    double var;
    double base = 2.0;
    int seed = -1;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -V FLOAT  variance [required]\n");
        fprintf(stderr, "         -s INT    rng seed [default: random]\n");
        fprintf(stderr, "         -C FLOAT  cover base [default: %.2f]\n", base);
        return -1;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);
    seed = read_int_arg(argc, argv, "-s", &seed);
    base = read_double_arg(argc, argv, "-C", &base);

    if (seed < 0)
    {
        random_device rd;
        default_random_engine gen(rd());
        uniform_int_distribution<int> dis{numeric_limits<int>::min(), numeric_limits<int>::max()};
        seed = dis(gen);
    }

    fprintf(stderr, "program: %s\ncommit: " GIT_COMMIT "\ndate and time: %s\n\n", program_str(argc, argv).c_str(), return_current_time_and_date().c_str());

    vector<Point> points = Point::random_points(n, var, seed);

    auto tree_start = mytimer::clock::now();

    SGTree tree(points, base);
    tree.build_tree();

    auto tree_end = mytimer::clock::now();
    double tree_time = mytimer::duration(tree_end - tree_start).count();

    fprintf(stderr, "[time=%.4f] :: (build_tree) [num_vertices=%lld,num_levels=%lld,base=%.2f]\n", tree_time, tree.num_vertices(), tree.num_levels(), base);

    tree.print_timing_results();

    return 0;
}

string program_str(int argc, char *argv[])
{
    stringstream ss;

    for (int i = 0; i < argc; ++i)
    {
        ss << argv[i] << " ";
    }

    return ss.str();
}

string return_current_time_and_date()
{
    /*
     * Shamelessly copied from here: https://stackoverflow.com/questions/17223096/outputting-date-and-time-in-c-using-stdchrono
     */

    auto now = chrono::system_clock::now();
    auto in_time_t = chrono::system_clock::to_time_t(now);

    stringstream ss;
    ss << put_time(localtime(&in_time_t), "%Y/%m/%d %X");
    return ss.str();
}
