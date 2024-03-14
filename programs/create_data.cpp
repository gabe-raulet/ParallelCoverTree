#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "VectorIO.h"
#include "read_args.h"

std::string get_filename(char *fname);

int main(int argc, char *argv[])
{
    size_t n;
    int d = 2;
    int seed = -1;
    int nthreads = 1;
    double var = 1.;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -o FILE   output filename [required]\n");
        fprintf(stderr, "         -d INT    point dimension [default: %d]\n", d);
        fprintf(stderr, "         -V FLOAT  rng variance [default: %.2f]\n", var);
        fprintf(stderr, "         -s INT    rng seed [default: random]\n");
        fprintf(stderr, "         -t INT    number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "         -h        help message\n");
        return 1;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    std::string fname = get_filename(read_string_arg(argc, argv, "-o", NULL));
    d = read_int_arg(argc, argv, "-d", &d);
    seed = read_int_arg(argc, argv, "-s", &seed);
    var = read_double_arg(argc, argv, "-V", &var);
    nthreads = std::min(read_int_arg(argc, argv, "-t", &nthreads), omp_get_max_threads());

    fprintf(stderr, "Parameters: [n=%lld, d=%d, seed=%d, filename='%s', variance=%.2f, nthreads=%d]\n", n, d, seed, fname.c_str(), var, nthreads);

    double t = -omp_get_wtime();

    std::random_device rd;
    std::vector<std::default_random_engine> gens;

    for (int i = 1; i <= nthreads; ++i)
    {
        gens.emplace_back(seed < 0? rd() : seed*17*i);
    }

    std::vector<float> pointmem(d*n);

    #pragma omp parallel firstprivate(d, n, var) num_threads(nthreads)
    {
        assert(nthreads == omp_get_num_threads());
        std::normal_distribution dis{0.0, std::sqrt(var)};
        auto& gen = gens[omp_get_thread_num()];

        #pragma omp for
        for (size_t i = 0; i < d*n; ++i)
            pointmem[i] = dis(gen);
    }

    t += omp_get_wtime();
    fprintf(stderr, "Generated %lld points of dimension %d in %.4f seconds\n", n, d, t);

    t = -omp_get_wtime();
    write_vecs_file(fname.c_str(), d, pointmem);
    t += omp_get_wtime();

    fprintf(stderr, "Wrote points to file '%s' in %.4f seconds\n", fname.c_str(), t);

    return 0;
}

std::string get_filename(char *fname)
{
    std::string filename(fname);
    char *dot = strrchr(fname, '.');

    if (dot && !strcmp(dot, ".fvecs"))
        return filename;

    filename += ".fvecs";
    return filename;
}
