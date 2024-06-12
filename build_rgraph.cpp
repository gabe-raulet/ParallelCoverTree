#include <vector>
#include <numeric>
#include <unistd.h>
#include "hrfsize.h"
#include "bforce.h"
#include "cover_tree.h"
#include "timers.h"
#include "ptraits.h"
#include "common.h"
#include "mpi_comm.h"

using namespace std;

using PointTraits = Traits<Real, POINT_DIM>;
using Point = PointTraits::Point;
using Distance = PointTraits::Distance;
using Index = int64_t;
using TBForce = BForce<PointTraits, double, Index>;
using TCoverTree = CoverTree<PointTraits, double, Index>;

using IndexSet = typename TBForce::IndexSet;
using IndexVector = typename TBForce::IndexVector;
using IndexSetVector = typename TBForce::IndexSetVector;

void read_options(int argc, char *argv[], char *&infname, char *&outprefix, double& cutoff, int& iters, double& damping, double& base, MPICommunicator& m_comm);
void read_points_file(vector<Point>& mypoints, const char *infname, MPICommunicator& m_comm);
void build_bforce_index(TBForce& bforce, const vector<Point>& mypoints, double cutoff);
void build_tree_index(TCoverTree& tree, const vector<Point>& mypoints, double cutoff, double base);
void build_bforce_rgraph(const TBForce& bforce, double radius, IndexSetVector& rgraph, int iter);
void build_tree_rgraph(const TCoverTree& tree, double radius, IndexSetVector& rgraph, int iter);
void write_rgraph_file(const IndexSetVector& rgraph, char *outprefix, int iter, MPI_Comm comm);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPICommunicator m_comm(MPI_COMM_WORLD);
    int myrank = m_comm.getrank();
    int nprocs = m_comm.getsize();

    char *infname, *outprefix = NULL;
    double cutoff;
    double damping = 0.5;
    double base = -1.0;
    int iters = 5;

    MPITimer timer(m_comm.getcomm());
    timer.start_timer();

    read_options(argc, argv, infname, outprefix, cutoff, iters, damping, base, m_comm);

    vector<Point> mypoints;
    vector<IndexSetVector> rgraphs;
    TBForce bforce(m_comm.getcomm());
    TCoverTree tree(m_comm.getcomm());

    read_points_file(mypoints, infname, m_comm);

    bool build_tree = (base > 0);

    if (build_tree) build_tree_index(tree, mypoints, cutoff, base);
    else build_bforce_index(bforce, mypoints, cutoff);

    double radius = cutoff;
    for (int i = 0; i < iters; ++i)
    {
        rgraphs.emplace_back();
        IndexSetVector& rgraph = rgraphs.back();
        if (build_tree) build_tree_rgraph(tree, radius, rgraph, i+1);
        else build_bforce_rgraph(bforce, radius, rgraph, i+1);
        radius *= damping;
    }

    if (outprefix)
        for (int iter = 0; const auto& rgraph : rgraphs)
            write_rgraph_file(rgraph, outprefix, ++iter, m_comm.getcomm());

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,nprocs=%d,msg::%s] command:", timer.maxtime, timer.avgtime, nprocs, __func__);
        for (int i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n");
    }

    MPI_Finalize();
    return 0;
}

void usage(const char *prg, int iters, double damping, MPI_Comm comm)
{
    int myrank;
    MPI_Comm_rank(comm, &myrank);

    if (!myrank)
    {
        fprintf(stderr, "Usage: %s [options] <ptsfname> <cutoff>\n", prg);
        fprintf(stderr, "Options: -n INT   number of iterations [%d]\n", iters);
        fprintf(stderr, "         -d FLOAT damping factor [%.2f]\n", damping);
        fprintf(stderr, "         -B FLOAT cover tree base (implies cover tree usage)\n");
        fprintf(stderr, "         -o FILE  output graph prefix\n");
        fprintf(stderr, "         -h       help message\n");
    }

    MPI_Finalize();
    exit(1);
}

void read_options(int argc, char *argv[], char *&infname, char *&outprefix, double& cutoff, int& iters, double& damping, double& base, MPICommunicator& m_comm)
{
    int myrank = m_comm.getrank();
    int nprocs = m_comm.getsize();

    int c;
    while ((c = getopt(argc, argv, "n:o:d:B:h")) >= 0)
    {
        if      (c == 'n') iters = atoi(optarg);
        else if (c == 'o') outprefix = optarg;
        else if (c == 'd') damping = atof(optarg);
        else if (c == 'B') base = atof(optarg);
        else if (c == 'h') usage(argv[0], iters, damping, m_comm.getcomm());
    }

    if (argc - optind < 2)
    {
        if (!myrank) fprintf(stderr, "[err::%s] missing argument(s)\n", __func__);
        usage(argv[0], iters, damping, m_comm.getcomm());
    }

    infname = argv[optind++];
    cutoff = atof(argv[optind]);
}

void read_points_file(vector<Point>& mypoints, const char *infname, MPICommunicator& m_comm)
{
    int myrank = m_comm.getrank();
    int nprocs = m_comm.getsize();

    MPITimer timer(m_comm.getcomm());
    timer.start_timer();

    vector<Point> points;
    vector<int> sendcounts;

    if (!myrank)
    {
        PointTraits::read_from_file(points, infname);
        sendcounts.resize(nprocs);
        fill(sendcounts.begin(), sendcounts.end(), points.size()/nprocs);
        sendcounts.back() = points.size() - (nprocs-1)*sendcounts.back();
    }

    m_comm.scatterv(points, sendcounts, mypoints, 0);

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] read %lu points from file '%s' [size=%s]\n", timer.maxtime, timer.avgtime, __func__, points.size(), infname, HumanReadable::str(infname).c_str());
    }
}

void build_bforce_index(TBForce& bforce, const vector<Point>& mypoints, double cutoff)
{
    int myrank, nprocs;
    MPI_Comm comm = bforce.getcomm();
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm);
    timer.start_timer();

    bforce.build(mypoints, cutoff);

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] brute force indexed %lld points with cutoff %.3f\n", timer.maxtime, timer.avgtime, __func__, bforce.gettotsize(), cutoff);
    }
}

void build_tree_index(TCoverTree& tree, const vector<Point>& mypoints, double cutoff, double base)
{
    int myrank, nprocs;
    MPI_Comm comm = tree.getcomm();
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm);
    timer.start_timer();

    tree.build(mypoints, cutoff, base);

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] cover tree indexed %lld points with cutoff %.3f\n", timer.maxtime, timer.avgtime, __func__, tree.gettotsize(), cutoff);
    }
}

void build_bforce_rgraph(const TBForce& bforce, double radius, IndexSetVector& rgraph, int iter)
{
    int myrank, nprocs;
    MPI_Comm comm = bforce.getcomm();
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm);
    timer.start_timer();

    Index my_n_edges = bforce.build_rgraph(radius, rgraph);

    timer.stop_timer();

    Index n_edges;
    MPI_Reduce(&my_n_edges, &n_edges, 1, MPI_INT64_T, MPI_SUM, 0, comm);

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] %d. built %.2f-graph [n_edges=%lld,avg_degree=%.2f]\n", timer.maxtime, timer.avgtime, __func__, iter, radius, n_edges, (n_edges+0.0)/bforce.gettotsize());
    }
}

void build_tree_rgraph(const TCoverTree& tree, double radius, IndexSetVector& rgraph, int iter)
{
    int myrank, nprocs;
    MPI_Comm comm = tree.getcomm();
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPITimer timer(comm);
    timer.start_timer();

    Index my_n_edges = tree.build_rgraph(radius, rgraph);

    timer.stop_timer();

    Index n_edges;
    MPI_Reduce(&my_n_edges, &n_edges, 1, MPI_INT64_T, MPI_SUM, 0, comm);

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] %d. built %.2f-graph [n_edges=%lld,avg_degree=%.2f]\n", timer.maxtime, timer.avgtime, __func__, iter, radius, n_edges, (n_edges+0.0)/tree.gettotsize());
    }
}

void write_rgraph_file(const IndexSetVector& rgraph, char *outprefix, int iter, MPI_Comm comm)
{
    MPICommunicator m_comm(comm);
    int myrank = m_comm.getrank();
    int nprocs = m_comm.getsize();

    MPITimer timer(comm);
    timer.start_timer();

    Index myoffset = 0;
    Index mysize = rgraph.size();
    m_comm.exscan(mysize, myoffset, MPI_SUM);

    stringstream ss;
    IndexVector neighs_v;
    Index my_n_edges = 0, n_edges;
    Index my_n_verts = rgraph.size(), n_verts;

    ss << outprefix << "_" << iter << ".rgraph";
    string fname = ss.str();

    ss.clear();
    ss.str({});

    for (Index u = 0; const auto& neighs : rgraph)
    {
        my_n_edges += neighs.size();
        neighs_v.assign(neighs.cbegin(), neighs.cend());

        for (Index v : neighs_v)
            ss << u+myoffset+1 << " " << v+1 << "\n";

        u++;
    }

    m_comm.reduce(my_n_edges, n_edges, 0, MPI_SUM);
    m_comm.reduce(my_n_verts, n_verts, 0, MPI_SUM);

    string sbuf;

    if (!myrank)
    {
        stringstream ss2;
        ss2 << n_verts << " " << n_edges << "\n" << ss.str();
        sbuf = ss2.str();
    }
    else sbuf = ss.str();

    vector<char> buf(sbuf.begin(), sbuf.end());

    MPI_Offset mybufsize = buf.size(), fileoffset = 0;
    m_comm.exscan(mybufsize, fileoffset, MPI_SUM);

    MPI_File fh;
    MPI_File_open(comm, fname.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at_all(fh, fileoffset, buf.data(), static_cast<int>(buf.size()), MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    timer.stop_timer();

    if (!myrank)
    {
        fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,msg::%s] wrote file '%s' [size=%s]\n", timer.maxtime, timer.avgtime, __func__, fname.c_str(), HumanReadable::str(fname.c_str()).c_str());
    }
}
