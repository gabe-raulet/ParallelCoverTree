#ifndef READ_ARGS_H_
#define READ_ARGS_H_

#include <stdint.h> /* int64_t */

int find_arg_idx(int argc, char *argv[], const char *argopt);
int read_int_arg(int argc, char *argv[], const char *argopt, int *defval);
int64_t read_formatted_int_arg(int argc, char *argv[], const char *argopt, int64_t *defval);
double read_double_arg(int argc, char  *argv[], const char *argopt, double *defval);
char* read_string_arg(int argc,  char *argv[], const char *argopt, char **defval);

#endif
