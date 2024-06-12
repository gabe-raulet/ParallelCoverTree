#include <iostream>
#include <random>
#include <vector>
#include <type_traits>
#include <unistd.h>
#include "hrfsize.h"
#include "timers.h"
#include "ptraits.h"
#include "common.h"

using namespace std;

using PointTraits = Traits<Real, POINT_DIM>;
using Point = PointTraits::Point;
using Distance = PointTraits::Distance;
using Index = int64_t;

void read_options(int argc, char *argv[], Index& n, char *&outfname, double& var, int& seed);
void generate_points(vector<Point>& points, double var, int seed);
void write_points_file(const vector<Point>& points, const char *outfname);

template <class Integer>
Integer readnum(char *str);

int main(int argc, char *argv[])
{
    Index n;
    char *outfname;
    double var = 10.0;
    int seed = -1;
    int c;

    LocalTimer timer;
    timer.start_timer();

    read_options(argc, argv, n, outfname, var, seed);

    vector<Point> points(n);

    generate_points(points, var, seed);
    write_points_file(points, outfname);

    timer.stop_timer();

    fprintf(stderr, "[time=%.3f,msg::%s] command:", timer.get_elapsed(), __func__);
    for (int i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");

    return 0;
}

void usage(const char *prg, double var, int seed)
{
    fprintf(stderr, "Usage: %s [options] <numpts> <outfname>\n", prg);
    fprintf(stderr, "Options: -V FLOAT variance [%.2f]\n", var);
    fprintf(stderr, "         -s INT   random seed [%s]\n", seed < 0? "random" : to_string(seed).c_str());
    fprintf(stderr, "         -h       help message\n");
    exit(1);
}

void read_options(int argc, char *argv[], Index& n, char *&outfname, double& var, int& seed)
{
    int c;
    while ((c = getopt(argc, argv, "V:s:h")) >= 0)
    {
        if      (c == 'V') var = atof(optarg);
        else if (c == 's') seed = atof(optarg);
        else if (c == 'h') usage(argv[0], var, seed);
    }

    if (argc - optind < 2)
    {
        fprintf(stderr, "[err::%s] missing argument(s)\n", __func__);
        usage(argv[0], var, seed);
    }

    n = readnum<Index>(argv[optind++]);
    outfname = argv[optind];
}

void generate_points(vector<Point>& points, double var, int seed)
{
    LocalTimer timer;
    timer.start_timer();

    random_device rd;
    if (seed < 0) seed = rd();
    default_random_engine gen(17*seed);
    normal_distribution<Real> dist{0, sqrt(static_cast<Real>(var))};
    PointTraits::fill_random_vec(points, gen, dist);

    timer.stop_timer();

    fprintf(stderr, "[time=%.3f,msg::%s] generated %lu random points [%s,var=%.2f,seed=%d]\n", timer.get_elapsed(), __func__, points.size(), PointTraits::name().c_str(), var, seed);
}

void write_points_file(const vector<Point>& points, const char *outfname)
{
    LocalTimer timer;
    timer.start_timer();

    PointTraits::write_to_file(points, outfname);

    timer.stop_timer();

    fprintf(stderr, "[time=%.3f,msg::%s] wrote %lu points to file '%s' [size=%s]\n", timer.get_elapsed(), __func__, points.size(), outfname, HumanReadable::str(outfname).c_str());
}

template <class Integer>
Integer readnum(char *str)
{
    double x;
    char *p;

    x = strtod(str, &p);

    if      (toupper(*p) == 'K') x *= (1LL << 10);
    else if (toupper(*p) == 'M') x *= (1LL << 20);
    else if (toupper(*p) == 'G') x *= (1LL << 30);

    return static_cast<Integer>(x + 0.499);
}
