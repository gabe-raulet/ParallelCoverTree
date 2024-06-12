#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include "hrfsize.h"
#include "timers.h"

using namespace std;

using Index = int64_t;
using IndexVector = vector<Index>;
using IndexVectorVector = vector<IndexVector>;

Index read_graph_file(IndexVectorVector& graph, const char *infname);
void compare_graphs(const IndexVectorVector& graph1, Index m1, const IndexVectorVector& graph2, Index m2);

int main(int argc, char *argv[])
{
    LocalTimer timer;
    timer.start_timer();

    if (argc != 3)
    {
        fprintf(stderr, "usage: %s <graph1> <graph2>\n", argv[0]);
        return 1;
    }

    IndexVectorVector graph1, graph2;
    Index m1, m2;

    m1 = read_graph_file(graph1, argv[1]);
    m2 = read_graph_file(graph2, argv[2]);

    compare_graphs(graph1, m1, graph2, m2);

    timer.stop_timer();

    fprintf(stderr, "[time=%.3f,msg::%s] command:", timer.get_elapsed(), __func__);
    for (int i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");

    return 0;
}

Index read_graph_file(IndexVectorVector& graph, const char *infname)
{
    LocalTimer timer;
    timer.start_timer();

    graph.clear();

    ifstream is;
    string line;
    Index n, m, u, v;

    is.open(infname);

    getline(is, line);
    istringstream header(line);
    header >> n >> m;

    graph.resize(n);

    while (getline(is, line))
    {
        istringstream ss(line);
        ss >> u >> v;
        graph[u-1].push_back(v-1);
    }

    is.close();

    timer.stop_timer();

    fprintf(stderr, "[time=%.3f,msg::%s] read file '%s' [num_vertices=%lld,num_edges=%lld,size=%s]\n", timer.get_elapsed(), __func__, infname, n, m, HumanReadable::str(infname).c_str());
    return m;
}

void compare_graphs(const IndexVectorVector& graph1, Index m1, const IndexVectorVector& graph2, Index m2)
{
    LocalTimer timer;
    timer.start_timer();

    Index n1 = graph1.size();
    Index n2 = graph2.size();

    if (n1 != n2 || m1 != m2)
    {
        timer.stop_timer();
        fprintf(stderr, "[time=%.3f,msg::%s] graphs differ\n", timer.get_elapsed(), __func__);
        return;
    }

    for (Index i = 0; i < n1; ++i)
        if (!is_permutation(graph1[i].cbegin(), graph1[i].cend(), graph2[i].cbegin()))
        {
            timer.stop_timer();
            fprintf(stderr, "[time=%.3f,msg::%s] graphs differ\n", timer.get_elapsed(), __func__);
            return;
        }

    timer.stop_timer();

    fprintf(stderr, "[time=%.3f,msg::%s] graphs are identical\n", timer.get_elapsed(), __func__);
}
