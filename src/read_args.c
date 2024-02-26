#include "read_args.h"
#include <stdio.h>
#include <string.h> /* strcmp */
#include <stdlib.h> /* atoi, atof, exit */
#include <ctype.h> /* toupper */

static inline int64_t read_formatted_int(char *str)
{
    double x;
    char *p;

    x = strtod(str, &p);

    if      (toupper(*p) == 'K') x *= (1LL<<10);
    else if (toupper(*p) == 'M') x *= (1LL<<20);
    else if (toupper(*p) == 'G') x *= (1LL<<30);

    return (int64_t)(x + 0.499);
}

int find_arg_idx(int argc, char *argv[], const char *argopt)
{
    for (int i = 1; i < argc; ++i)    
        if (!strcmp(argv[i], argopt))
            return i;

    return -1;
}

int read_int_arg(int argc, char *argv[], const char *argopt, int *defval)
{
    int arg_idx = find_arg_idx(argc, argv, argopt);

    if (arg_idx >= 0 && arg_idx+1 < argc)
        return atoi(argv[arg_idx+1]);

    if (!defval)
    {
        fprintf(stderr, "error: missing required argment '%s'\n", argopt);
        exit(-1);
    }

    return *defval;
}

int64_t read_formatted_int_arg(int argc, char *argv[], const char *argopt, int64_t *defval)
{
    int arg_idx = find_arg_idx(argc, argv, argopt);

    if (arg_idx >= 0 && arg_idx+1 < argc)
        return read_formatted_int(argv[arg_idx+1]);

    if (!defval)
    {
        fprintf(stderr, "error: missing required argment '%s'\n", argopt);
        exit(-1);
    }

    return *defval;
}

double read_double_arg(int argc, char *argv[], const char *argopt, double *defval)
{
    int arg_idx = find_arg_idx(argc, argv, argopt);

    if (arg_idx >= 0 && arg_idx+1 < argc)
        return atof(argv[arg_idx+1]);

    if (!defval)
    {
        fprintf(stderr, "error: missing required argment '%s'\n", argopt);
        exit(-1);
    }

    return *defval;
}

char* read_string_arg(int argc, char *argv[], const char *argopt, char **defval)
{
    int arg_idx = find_arg_idx(argc, argv, argopt);

    if (arg_idx >= 0 && arg_idx+1 < argc)
        return argv[arg_idx+1];


    if (!defval)
    {
        fprintf(stderr, "error: missing required argment '%s'\n", argopt);
        exit(-1);
    }

    return *defval;
}
