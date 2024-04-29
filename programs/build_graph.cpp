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

void write_graph_file(const vector<vector<int64_t>>& graph, const char *fname, int64_t n_edges = -1);

int main(int argc, char *argv[])
{
    double radius;
    char *ifname = NULL;
    char *ofname = NULL;
    double base = 2.0;
    bool verbose = false;
    bool skip_graph = false;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -i FILE   input filename [required]\n");
        fprintf(stderr, "         -r FLOAT  radius [required]\n");
        fprintf(stderr, "         -C FLOAT  cover base [default: %.2f]\n", base);
        fprintf(stderr, "         -o FILE   output filename\n");
        fprintf(stderr, "         -S        skip graph construction\n");
        fprintf(stderr, "         -v        verbose\n");
        fprintf(stderr, "         -h        help message\n");
        return -1;
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

    fprintf(stderr, "program: %s\ncommit: " GIT_COMMIT "\ndate and time: %s\n\n", program_str(argc, argv).c_str(), return_current_time_and_date().c_str());

    auto read_file_start = mytimer::clock::now();

    vector<Point> points = Point::from_file(ifname);

    auto read_file_end = mytimer::clock::now();
    double read_file_time = mytimer::duration(read_file_end - read_file_start).count();

    fprintf(stderr, "[time=%.4f] :: (read_points) [n=%lu,filename='%s']\n", read_file_time, points.size(), ifname);

    auto tree_start = mytimer::clock::now();

    CoverTree tree(points, base);
    tree.build_tree(verbose);

    auto tree_end = mytimer::clock::now();
    double tree_time = mytimer::duration(tree_end - tree_start).count();

    tree.print_timing_results();
    fprintf(stderr, "[time=%.4f] :: (build_tree) [num_vertices=%lld,num_levels=%lld,base=%.2f]\n", tree_time, tree.num_vertices(), tree.num_levels(), base);

    if (skip_graph) return 0;

    auto graph_start = mytimer::clock::now();

    auto graph = tree.build_epsilon_graph(radius);

    auto graph_end = mytimer::clock::now();
    double graph_time = mytimer::duration(graph_end - graph_start).count();

    int64_t n_edges = 0;
    for_each(graph.begin(), graph.end(), [&](const auto& item) { n_edges += item.size(); });

    fprintf(stderr, "[time=%.4f] :: (build_graph) [num_vertices=%lu,num_edges=%lld,avg_deg=%.4f,radius=%.3f]\n", graph_time, graph.size(), n_edges, (n_edges+0.0)/graph.size(), radius);
    fprintf(stderr, "[time=%.4f] :: (cover_tree_graph)\n", graph_time + tree_time);

    if (ofname)
    {
        auto write_graph_start = mytimer::clock::now();

        write_graph_file(graph, ofname, n_edges);

        auto write_graph_end = mytimer::clock::now();
        double write_graph_time = mytimer::duration(write_graph_end - write_graph_start).count();

        fprintf(stderr, "[time=%.4f] :: (write_graph) [filename='%s']\n", write_graph_time, ofname);
    }

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

void write_graph_file(const vector<vector<int64_t>>& graph, const char *fname, int64_t n_edges)
{
    int64_t n_verts = graph.size();
    if (n_edges < 0) for_each(graph.begin(), graph.end(), [&](const auto& item) { n_edges += item.size(); });

    FILE *f = fopen(fname, "w");
    fprintf(f, "%lld %lld\n", n_verts, n_edges);

    for (int64_t u = 0; u < n_verts; ++u)
    {
        vector<int64_t> dests = graph[u];
        sort(dests.begin(), dests.end());

        for (int64_t v : dests)
        {
            fprintf(f, "%lld %lld\n", u+1, v+1);
        }
    }

    fclose(f);
}
