#include "Point.h"
#include "CoverTree.h"
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
    double radius;
    char *fname = NULL;
    double base = 2.0;
    bool verbose = false;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -i FILE   input filename [required]\n");
        fprintf(stderr, "         -r FLOAT  radius [required]\n");
        fprintf(stderr, "         -C FLOAT  cover base [default: %.2f]\n", base);
        fprintf(stderr, "         -v        verbose\n");
        fprintf(stderr, "         -h        help message\n");
        return -1;
    }

    radius = read_double_arg(argc, argv, "-r", NULL);
    fname = read_string_arg(argc, argv, "-i", NULL);
    base = read_double_arg(argc, argv, "-C", &base);
    verbose = (find_arg_idx(argc, argv, "-v") >= 0);

    fprintf(stderr, "program: %s\ncommit: " GIT_COMMIT "\ndate and time: %s\n\n", program_str(argc, argv).c_str(), return_current_time_and_date().c_str());

    auto read_file_start = mytimer::clock::now();

    vector<Point> points = Point::from_file(fname);

    auto read_file_end = mytimer::clock::now();
    double read_file_time = mytimer::duration(read_file_end - read_file_start).count();

    int64_t n = points.size();
    fprintf(stderr, "[time=%.4f] :: (read_file) [n=%lld,filename='%s']\n", read_file_time, n, fname);

    Point::to_file(points, "test");

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
