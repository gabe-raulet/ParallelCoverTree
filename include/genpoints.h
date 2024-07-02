#ifndef GEN_POINTS_H_
#define GEN_POINTS_H_

#include <vector>
#include <random>
#include <algorithm>
#include "mpi_env.h"
#include "misc.h"

using namespace std;

template <class PointTraits>
class BlockPointGenerator
{
    public:

        using Point = typename PointTraits::Point;
        using Real = typename PointTraits::Real;
        using Comm = MPIEnv::Comm;

        BlockPointGenerator(int seed, int blocks);
        BlockPointGenerator(int seed, int blocks, Comm comm);

        void generate_points(vector<Point>& points, size_t totsize, double var);

    private:

        int blocks;
        Comm comm;
        vector<int> seeds;

        template <class Iter>
        void generate_block(Iter first, Iter last, int seed, double var);
};

#include "genpoints.hpp"

#endif
