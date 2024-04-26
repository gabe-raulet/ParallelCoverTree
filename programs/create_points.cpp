#include "Point.h"
#include "MyTimer.h"
#include "read_args.h"
#include "version.h"
#include <iostream>
#include <sstream>
#include <chrono>
#include <assert.h>
#include <stdio.h>

using namespace std;

string return_current_time_and_date();
string program_str(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    int64_t n;
    double var = 10.0;
    int seed = -1;
    char *fname = NULL;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -o FILE   output filename [required]\n");
        fprintf(stderr, "         -V FLOAT  variance [default: %.2f]\n", var);
        fprintf(stderr, "         -s INT    rng seed [default: random]\n");
        fprintf(stderr, "         -h        help message\n");
        return -1;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    var = read_double_arg(argc, argv, "-V", &var);
    seed = read_int_arg(argc, argv, "-s", &seed);
    fname = read_string_arg(argc, argv, "-o", NULL);

    if (seed < 0)
    {
        random_device rd;
        default_random_engine gen(rd());
        uniform_int_distribution<int> dis{0, numeric_limits<int>::max()};
        seed = dis(gen);
    }

    fprintf(stderr, "program: %s\ncommit: " GIT_COMMIT "\ndate and time: %s\n\n", program_str(argc, argv).c_str(), return_current_time_and_date().c_str());


    /*
     * Construct random point set
     */
    auto random_start = mytimer::clock::now();

    vector<Point> points = Point::random_points(n, var, seed);

    auto random_end = mytimer::clock::now();
    double random_time = mytimer::duration(random_end - random_start).count();

    fprintf(stderr, "[time=%.4f] :: (random_points) [n=%lld,var=%.2f,seed=%d]\n", random_time, n, var, seed);

    auto io_start = mytimer::clock::now();

    Point::to_file(points, fname);

    auto io_end = mytimer::clock::now();
    double io_time = mytimer::duration(io_end - io_start).count();

    fprintf(stderr, "[time=%.4f] :: (to_file) [filename='%s']\n", io_time, fname);

    return 0;
}

string program_str(int argc, char *argv[])
{
    stringstream ss;

    for (int i = 0; i < argc; ++i)
        ss << argv[i] << " ";

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
