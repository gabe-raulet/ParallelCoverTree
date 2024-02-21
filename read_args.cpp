#include "read_args.h"
#include <string>
#include <cstring>

int find_int_arg(int argc, char *argv[], const char *option, int default_value)
{
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc-1)
        return std::stoi(argv[iplace+1]);

    return default_value;
}

int find_arg_idx(int argc, char *argv[], const char *option)
{
    for (int i = 1; i < argc; ++i)
        if (!strcmp(argv[i], option))
            return i;

    return -1;
}

char* find_string_option(int argc, char *argv[], const char *option, char *default_value)
{
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc-1)
        return argv[iplace+1];

    return default_value;
}
