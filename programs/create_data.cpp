#include <algorithm>
#include <vector>
#include <iostream>
#include <random>
#include <tuple>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "read_args.h"

int main(int argc, char *argv[])
{
    size_t n;
    int d = 2;
    int seed = -1;
    double minval = -1.0;
    double maxval = 1.0;
    int nthreads = 1;
    char *filename = NULL;

    if (argc == 1 || find_arg_idx(argc, argv, "-h") >= 0)
    {
        fprintf(stderr, "Usage: %s [options]\n", argv[0]);
        fprintf(stderr, "Options: -n INT    number of points [required]\n");
        fprintf(stderr, "         -o FILE   output filename [required]\n");
        fprintf(stderr, "         -d INT    point dimension [default: %d]\n", d);
        fprintf(stderr, "         -L FLOAT  min element value [default: %.2f]\n", minval);
        fprintf(stderr, "         -U FLOAT  max element value [default: %.2f]\n", maxval);
        fprintf(stderr, "         -s INT    rng seed [default: random]\n");
        fprintf(stderr, "         -t INT    number of threads [default: %d]\n", nthreads);
        fprintf(stderr, "         -h        help message\n");
        return 1;
    }

    n = read_formatted_int_arg(argc, argv, "-n", NULL);
    d = read_int_arg(argc, argv, "-d", &d);
    minval = read_double_arg(argc, argv, "-L", &minval);
    maxval = read_double_arg(argc, argv, "-U", &maxval);
    seed = read_int_arg(argc, argv, "-s", &seed);
    filename = read_string_arg(argc, argv, "-o", NULL);
    nthreads = std::min(read_int_arg(argc, argv, "-t", &nthreads), omp_get_max_threads());

    assert(minval <= maxval);

    double t = -omp_get_wtime();

    std::random_device rd;
    std::vector<std::default_random_engine> gens;

    for (int i = 1; i <= nthreads; ++i)
    {
        gens.emplace_back(seed < 0? rd() : seed*17*i);
    }

    std::vector<float> pointmem(d*n);

    #pragma omp parallel firstprivate(d, n, minval, maxval) num_threads(nthreads)
    {
        assert(nthreads == omp_get_num_threads());
        std::uniform_real_distribution dis(minval, maxval);
        auto& gen = gens[omp_get_thread_num()];

        #pragma omp for
        for (size_t i = 0; i < d*n; ++i)
            pointmem[i] = dis(gen);
    }

    const float *p = pointmem.data();
    FILE *f = fopen(filename, "w");

    for (size_t i = 0; i < n; ++i)
    {
        fwrite(&d, sizeof(int), 1, f);
        fwrite(&p[i*d], sizeof(float), d, f);
    }

    fclose(f);

    t += omp_get_wtime();
    fprintf(stderr, "Finished writing %lld points of dimension %d in %.4f seconds\n", n, d, t);

    return 0;
}
