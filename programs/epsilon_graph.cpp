#include "Point.h"
#include "CoverTree.h"
#include "read_args.h"
#include <iostream>
#include <chrono>
#include <assert.h>
#include <stdio.h>

using namespace std;

template <class T>
struct MyTimer
{
    using clock = chrono::high_resolution_clock;
    using duration = chrono::duration<T>;
};

using mytimer = MyTimer<double>;

bool graphs_equal(const vector<vector<int64_t>> g1, const vector<vector<int64_t>> g2);

int main(int argc, char *argv[])
{
    int64_t n;
    double var;
    double radius;
    double base = 2.0;
    int seed = -1;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -V FLOAT  variance [required]\n");
        fprintf(stderr, "         -r FLOAT  radius [required]\n");
        fprintf(stderr, "         -s INT    rng seed [default: random]\n");
        fprintf(stderr, "         -C FLOAT  cover base [default: %.2f]\n", base);
        return -1;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    var = read_double_arg(argc, argv, "-V", NULL);
    radius = read_double_arg(argc, argv, "-r", NULL);
    seed = read_int_arg(argc, argv, "-s", &seed);
    base = read_double_arg(argc, argv, "-C", &base);

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

    auto tree_end = mytimer::clock::now();
    double tree_time = mytimer::duration(tree_end - tree_start).count();

    fprintf(stderr, "[time=%.4f] :: (build_tree) [num_vertices=%lld,num_levels=%lld,base=%.2f]\n", tree_time, tree.num_vertices(), tree.num_levels(), base);

    /*
     * Build epsilon graph from cover tree
     */

    auto graph_start = mytimer::clock::now();

    auto epsilon_graph = tree.build_epsilon_graph(radius);

    auto graph_end = mytimer::clock::now();
    double graph_time = mytimer::duration(graph_end - graph_start).count();

    int64_t n_edges = 0;
    for_each(epsilon_graph.begin(), epsilon_graph.end(), [&](const auto& item) { n_edges += item.size(); });

    fprintf(stderr, "[time=%.4f] :: (build_graph) [num_vertices=%lu,num_edges=%lld,avg_deg=%.4f]\n", graph_time, epsilon_graph.size(), n_edges, (n_edges+0.0)/epsilon_graph.size());
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

