#include "VectorIO.h"
#include <filesystem>
#include <iostream>
#include <stdio.h>

namespace fs = std::filesystem;

void write_vecs_file(const char *fname, int dim, const std::vector<float>& pdata)
{
    size_t n = pdata.size() / dim;
    const float *p = pdata.data();
    FILE *f;

    f = fopen(fname, "wb");

    for (size_t i = 0; i < n; ++i)
    {
        fwrite(&dim, sizeof(int), 1, f);
        fwrite(&p[i*dim], sizeof(float), dim, f);
    }

    fclose(f);
}

std::vector<float> read_vecs_file(const char *fname, int *dim, size_t *npts)
{
    std::vector<float> pdata;
    size_t filesize, n;
    int d;
    FILE *f;
    fs::path path = fname;

    if (!fs::exists(path) || !fs::is_regular_file(path))
    {
        std::cerr << "error: unable to open " << std::quoted(fname) << std::endl;
        exit(1);
    }

    filesize = fs::file_size(path);

    f = fopen(fname, "rb");
    fread(&d, sizeof(int), 1, f);
    n = filesize / (4*(d+1));
    fseek(f, SEEK_SET, 0);

    if (n == 0 || d <= 0)
    {
        std::cerr << "error: unable to open " << std::quoted(fname) << std::endl;
        exit(1);
    }

    pdata.resize(n*d);
    float *p = pdata.data();

    for (size_t i = 0; i < n; ++i)
    {
        fread(&d, sizeof(int), 1, f);
        fread(&p[i*d], sizeof(float), d, f);
    }

    fclose(f);

    if (dim) *dim = d;
    if (npts) *npts = n;

    return pdata;
}
