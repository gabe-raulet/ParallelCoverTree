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

bool graphs_equal(const vector<vector<int64_t>> g1, const vector<vector<int64_t>> g2);

string return_current_time_and_date();
string program_str(int argc, char *argv[]);

int main(int argc, char *argv[])
{
    int64_t n;
    double var;
    double radius;
    double base = 2.0;
    int seed = -1;
    int skip_test = 0;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -V FLOAT  variance [required]\n");
        fprintf(stderr, "         -r FLOAT  radius [required]\n");
        fprintf(stderr, "         -s INT    rng seed [default: random]\n");
        fprintf(stderr, "         -C FLOAT  cover base [default: %.2f]\n", base);
        fprintf(stderr, "         -F        skip test\n");
        fprintf(stderr, "         -h        help message\n");
        return -1;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);
    radius = read_double_arg(argc, argv, "-r", NULL);
    seed = read_int_arg(argc, argv, "-s", &seed);
    base = read_double_arg(argc, argv, "-C", &base);
    skip_test = !!(find_arg_idx(argc, argv, "-F") >= 0);

    if (seed < 0)
    {
        random_device rd;
        default_random_engine gen(rd());
        uniform_int_distribution<int> dis{numeric_limits<int>::min(), numeric_limits<int>::max()};
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

    /*
     * Build cover tree
     */

    auto tree_start = mytimer::clock::now();

    CoverTree tree(points, base);
    tree.build_tree();

    auto tree_end = mytimer::clock::now();
    double tree_time = mytimer::duration(tree_end - tree_start).count();

    tree.print_timing_results();
    fprintf(stderr, "[time=%.4f] :: (build_tree) [num_vertices=%lld,num_levels=%lld,base=%.2f]\n", tree_time, tree.num_vertices(), tree.num_levels(), base);

    if (skip_test) return 0;

    /*
     * Build epsilon graph from cover tree
     */

    auto graph_start = mytimer::clock::now();

    auto epsilon_graph = tree.build_epsilon_graph(radius);

    auto graph_end = mytimer::clock::now();
    double graph_time = mytimer::duration(graph_end - graph_start).count();

    int64_t n_edges = 0;
    for_each(epsilon_graph.begin(), epsilon_graph.end(), [&](const auto& item) { n_edges += item.size(); });

    fprintf(stderr, "[time=%.4f] :: (build_graph) [num_vertices=%lu,num_edges=%lld,avg_deg=%.4f,radius=%.3f]\n", graph_time, epsilon_graph.size(), n_edges, (n_edges+0.0)/epsilon_graph.size(), radius);
    fprintf(stderr, "[time=%.4f] :: (cover_tree_graph)\n", graph_time + tree_time);

    /*
     * Brute force construction of epsilon graph
     */

    auto bf_start = mytimer::clock::now();

    vector<vector<int64_t>> graph(points.size());

    for (int64_t u = 0; u < points.size(); ++u)
        for (int64_t v = 0; v < points.size(); ++v)
            if (points[u].distance(points[v]) <= radius)
                graph[u].push_back(v);

    auto bf_end = mytimer::clock::now();
    double bf_time = mytimer::duration(bf_end - bf_start).count();

    int64_t m = 0;
    for_each(graph.begin(), graph.end(), [&](const auto& item) { m += item.size(); });

    fprintf(stderr, "[time=%.4f] :: (brute_force_graph)\n", bf_time);
    fprintf(stderr, "\n");

    /*
     * Check that cover tree graph equals brute force graph
     */

    if (graphs_equal(graph, epsilon_graph))
    {
        fprintf(stderr, "Graph construction was successful\n");
    }
    else
    {
        fprintf(stderr, "Graph construction was NOT successful\n");
    }

    double speedup = (bf_time + 1e-16) / (graph_time + tree_time + 1e-16);
    fprintf(stderr, "Cover tree vs. brute force graph construction: %.2fx %s\n", speedup <= 1.? 1./speedup : speedup, speedup <= 1.? "slowdown" : "speedup");

    return 0;
}

bool graphs_equal(const vector<vector<int64_t>> g1, const vector<vector<int64_t>> g2)
{
    if (g1.size() != g2.size()) return false;

    int64_t n = g1.size();
    for (int64_t u = 0; u < n; ++u)
    {
        vector<int64_t> r1 = g1[u];
        vector<int64_t> r2 = g2[u];
        sort(r1.begin(), r1.end());
        sort(r2.begin(), r2.end());
        if (r1 != r2) return false;
    }

    return true;
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
