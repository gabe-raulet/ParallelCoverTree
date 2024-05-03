#include "Point.h"
#include "DistCoverTree.h"
#include "MPITimer.h"
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
void dist_write_graph_file(const vector<vector<int64_t>>& my_edges, const char *fname, MPI_Comm comm);

int main(int argc, char *argv[])
{
    int myrank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    double radius;
    char *ifname = NULL;
    char *ofname = NULL;
    double base = 2.0;
    bool verbose = false;
    bool skip_graph = false;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        if (!myrank)
        {
            fprintf(stderr, "Usage: %s [options]\n", argv[0]);
            fprintf(stderr, "Options: -i FILE   input filename [required]\n");
            fprintf(stderr, "         -r FLOAT  radius [required]\n");
            fprintf(stderr, "         -C FLOAT  cover base [default: %.2f]\n", base);
            fprintf(stderr, "         -o FILE   output filename\n");
            fprintf(stderr, "         -S        skip graph construction\n");
            fprintf(stderr, "         -v        verbose\n");
            fprintf(stderr, "         -h        help message\n");
        }

        MPI_Finalize();
        return 0;
    }

    ifname = read_string_arg(argc, argv, "-i", NULL);
    base = read_double_arg(argc, argv, "-C", &base);
    verbose = (find_arg_idx(argc, argv, "-v") >= 0);
    skip_graph = (find_arg_idx(argc, argv, "-S") >= 0);

    if (!skip_graph)
    {
        radius = read_double_arg(argc, argv, "-r", NULL);
    }

    if (find_arg_idx(argc, argv, "-o") >= 0)
    {
        ofname = read_string_arg(argc, argv, "-o", NULL);
    }

    if (!myrank) fprintf(stderr, "program: %s\ncommit: " GIT_COMMIT "\ndate and time: %s\nMPI ranks: %d\n\n", program_str(argc, argv).c_str(), return_current_time_and_date().c_str(), nprocs);

    MPITimer timer(MPI_COMM_WORLD, 0);

    timer.start_timer();

    vector<Point> points, mypoints;

    if (!myrank) points = Point::from_file(ifname);
    mypoints = Point::scatter(points, 0, MPI_COMM_WORLD);

    timer.stop_timer();

    if (!myrank) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (read_points) [n=%lu,filename='%s']\n", timer.get_max_time(), timer.get_avg_time(), points.size(), ifname);

    timer.start_timer();

    DistCoverTree tree(mypoints, base, MPI_COMM_WORLD);
    tree.build_tree(verbose);

    timer.stop_timer();
    tree.print_timing_results();

    double graph_time = timer.get_max_time();

    if (!myrank) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (dist_build_tree) [num_vertices=%lld,num_levels=%lld,base=%.2f]\n", timer.get_max_time(), timer.get_avg_time(), tree.num_vertices(), tree.num_levels(), base);

    if (skip_graph)
    {
        MPI_Finalize();
        return 0;
    }

    timer.start_timer();

    vector<vector<int64_t>> my_edges = tree.build_epsilon_graph(radius);

    timer.stop_timer();
    graph_time += timer.get_max_time();

    int64_t my_n_edges = 0, n_edges;
    for_each(my_edges.begin(), my_edges.end(), [&](const auto& item) { my_n_edges += item.size(); });
    MPI_Reduce(&my_n_edges, &n_edges, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!myrank) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (dist_build_graph) [num_vertices=%lu,num_edges=%lld,avg_deg=%.4f,radius=%.03f]\n", timer.get_max_time(), timer.get_avg_time(), points.size(), n_edges, (n_edges+0.0)/points.size(), radius);
    if (!myrank) fprintf(stderr, "[tottime=%.4f] :: (dist_cover_tree_graph)\n", graph_time);

    if (ofname)
    {
        timer.start_timer();
        dist_write_graph_file(my_edges, ofname, MPI_COMM_WORLD);
        timer.stop_timer();

        if (!myrank) fprintf(stderr, "[maxtime=%.4f,avgtime=%.4f] :: (dist_write_graph) [filename='%s']\n", timer.get_max_time(), timer.get_avg_time(), ofname);
    }

    MPI_Finalize();
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

void dist_write_graph_file(const vector<vector<int64_t>>& my_edges, const char *fname, MPI_Comm comm)
{
    int myrank, nprocs;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    int64_t my_n_edges = 0;
    int64_t my_n_verts = my_edges.size();
    int64_t myoffset;

    for_each(my_edges.cbegin(), my_edges.cend(), [&](const auto& item) { my_n_edges += item.size(); });

    MPI_Exscan(&my_n_verts, &myoffset, 1, MPI_INT64_T, MPI_SUM, comm);
    if (!myrank) myoffset = 0;

    int64_t n_edges, n_verts;
    MPI_Reduce(&my_n_verts, &n_verts, 1, MPI_INT64_T, MPI_SUM, 0, comm);
    MPI_Reduce(&my_n_edges, &n_edges, 1, MPI_INT64_T, MPI_SUM, 0, comm);

    stringstream ss;

    if (!myrank)
    {
        ss << n_verts << " " << n_edges << "\n";
    }

    for (int64_t u = 0; u < my_n_verts; ++u)
    {
        vector<int64_t> dests = my_edges[u];
        sort(dests.begin(), dests.end());

        for (int64_t v : dests)
        {
            ss << (u+myoffset+1) << " " << (v+1) << "\n";
        }
    }

    string sbuf = ss.str();
    vector<char> buf(sbuf.begin(), sbuf.end());

    MPI_Offset mysize = buf.size(), fileoffset;
    MPI_Exscan(&mysize, &fileoffset, 1, MPI_OFFSET, MPI_SUM, comm);
    if (!myrank) fileoffset = 0;

    MPI_File fh;
    MPI_File_open(comm, fname, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at_all(fh, fileoffset, buf.data(), static_cast<int>(buf.size()), MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}
