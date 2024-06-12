#ifndef PTRAITS_H_
#define PTRAITS_H_

#include <array>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <type_traits>

using namespace std;

/*
 * The Traits class defines all the necessary information and routines
 * to work with point data that can be represented as fixed-size array
 * of floating-point numbers. The template parameters Real and D define
 * the underlying floating point type (either float or double) and the
 * dimensionality of points, respectively.
 *
 * Particularly important are the distance() function, which returns
 * a Distance functor that can be used to compute the metric distance
 * between two points. A type declaration for the Point type is also
 * provided.
 *
 * Also useful are the fill_random routines (for generating random points
 * with user-provided RNG distributions and generators, as well as
 * methods for reading and writing points to disk.
 */

template <class Real, int D>
struct Traits
{
    static_assert(is_same_v<Real, float> || is_same_v<Real, double>);

    using Point = array<Real, D>;

    struct Distance { Real operator()(const Point& p, const Point& q); };
    static Distance distance() { return Distance(); }
    static int dimension() { return D; }
    static string name();

    template <class RandomGen, class RandomDist>
    static void fill_random(Point& point, RandomGen& gen, RandomDist& dist);

    template <class RandomGen, class RandomDist>
    static void fill_random_vec(vector<Point>& points, RandomGen& gen, RandomDist& dist);

    static void write_to_file(const vector<Point>& points, const char *outfname);
    static void read_from_file(vector<Point>& points, const char *infname);
};

#include "ptraits.hpp"

#endif
