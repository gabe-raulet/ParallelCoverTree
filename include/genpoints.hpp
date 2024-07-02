template <class PointTraits>
BlockPointGenerator<PointTraits>::BlockPointGenerator(int seed, int blocks) :
    blocks(blocks),
    comm(Comm::comm_self())
{
    random_device rd;
    if (seed < 0) seed = rd();
    else seed = ((seed+0xde)*0xadb)-0xeef;

    seeds.resize(blocks);

    for (int i = 0; i < blocks; ++i)
    {
        seeds[i] = 17*i*seed;
    }
}

template <class PointTraits>
BlockPointGenerator<PointTraits>::BlockPointGenerator(int seed, int blocks, Comm comm) :
    BlockPointGenerator(seed, blocks)
{
    this->comm = comm;
    assert(blocks >= comm.size());
    comm.bcast(seeds, 0);
}

template <class PointTraits>
template <class Iter>
void BlockPointGenerator<PointTraits>::generate_block(Iter first, Iter last, int seed, double var)
{
    default_random_engine gen(seed);
    normal_distribution<Real> dist{0.0, sqrt(static_cast<Real>(var))};
    PointTraits::fill_random_points(first, last, gen, dist);
}

template <class PointTraits>
void BlockPointGenerator<PointTraits>::generate_points(vector<Point>& points, size_t totsize, double var)
{
    int myrank = comm.rank();
    int nprocs = comm.size();

    vector<int> ptcnts(blocks); /* points per block */
    vector<int> blkcnts(nprocs), blkdisps(nprocs); /* blocks per rank and their displacements */
    vector<int>  ptsprank(nprocs, 0); /* points per rank */

    get_balanced_counts(ptcnts, totsize);
    get_balanced_counts(blkcnts, blocks);

    reverse(blkcnts.begin(), blkcnts.end());

    int i, j, k = 0;

    for (i = 0; i < nprocs; ++i)
        for (j = 0; j < blkcnts[i]; ++j)
            ptsprank[i] += ptcnts[k++];

    points.resize(ptsprank[myrank]);
    exclusive_scan(blkcnts.begin(), blkcnts.end(), blkdisps.begin(), static_cast<int>(0));

    int first = blkdisps[myrank];
    int last = first + blkcnts[myrank];
    auto it = points.begin();

    for (i = first; i < last; ++i)
    {
        generate_block(it, it + ptcnts[i], seeds[i], var);
        it += ptcnts[i];
    }
}
