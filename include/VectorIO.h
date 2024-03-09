#ifndef VECTOR_IO_H_
#define VECTOR_IO_H_

#include <vector>
#include <stddef.h>

std::vector<float> read_vecs_file(const char *fname, int *dim, size_t *npts);
void write_vecs_file(const char *fname, int dim, const std::vector<float>& pdata);

#endif
