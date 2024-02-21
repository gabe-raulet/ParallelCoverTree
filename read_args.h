#ifndef READ_ARGS_H_
#define READ_ARGS_H_

int find_arg_idx(int argc, char *argv[], const char *option);
int find_int_arg(int argc, char *argv[], const char *option, int default_value);
char* find_string_option(int argc, char *argv[], const char *option, char *default_value);

#endif
