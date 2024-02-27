#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
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

    fprintf(stderr, "         -n INT    number of points [%lld]\n", n);
    fprintf(stderr, "         -d INT    point dimension [%d]\n", d);
    fprintf(stderr, "         -r FLOAT  start radius [%.4f]\n", r0);
    fprintf(stderr, "         -i FLOAT  radius increment [%.4f]\n", inc);
    fprintf(stderr, "         -f INT    filter size [%d]\n", filtsize);
    fprintf(stderr, "         -s INT    seed [%d]\n", seed);
    fprintf(stderr, "         -o STR    output prefix ['%s']\n", output);
    fprintf(stderr, "         -h        help message\n");

    return 0;
}
