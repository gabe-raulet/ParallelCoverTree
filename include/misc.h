#ifndef MISC_H_
#define MISC_H_

#include <filesystem>
#include <sstream>
#include <fstream>
#include <string>
#include <concepts>
#include <stdio.h>

using namespace std;

size_t get_file_size(const char *fname)
{
    ifstream is;
    size_t fsize;

    is.open(fname, ios::binary | ios::in);
    is.seekg(0, is.end);
    fsize = is.tellg();
    is.seekg(0, is.beg);

    is.close();
    return fsize;
}

struct PrettyFileSize
{
    uintmax_t size;

    PrettyFileSize(const char *fname) : size(get_file_size(fname)) {}

    static string str(const char *fname)
    {
        stringstream ss;
        ss << PrettyFileSize(fname);
        return ss.str();
    }

    template <class c, class t>
    friend basic_ostream<c,t>& operator<<(basic_ostream<c,t>& os, PrettyFileSize hr)
    {
        int i{};
        double mantissa = hr.size;
        for (; mantissa >= 1024.0; mantissa /= 1024.0, ++i);
        os << ceil(mantissa * 10.) / 10. << i["BKMGTPE"];
        return i? os << "B" : os;
    }
};

void main_msg(int argc, char *argv[], double elapsed)
{
    fprintf(stderr, "[time=%.3f,msg::main] command:", elapsed);
    for (int i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");
}

void main_msg(int argc, char *argv[], double maxtime, double sumtime, int nprocs)
{
    fprintf(stderr, "[maxtime=%.3f,avgtime=%.3f,nprocs=%d,msg::main] command:", maxtime, sumtime/nprocs, nprocs);
    for (int i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");
}

template <integral T>
void get_balanced_counts(vector<T>& counts, size_t totsize)
{
    T blocks = counts.size();
    fill(counts.begin(), counts.end(), totsize/blocks);

    counts.back() = totsize - (blocks-1)*(totsize/blocks);

    T diff = counts.back() - counts.front();

    for (T i = 0; i < diff; ++i)
    {
        counts[blocks-1-i]++;
        counts[blocks-1]--;
    }
}

template <integral T>
T read_integer(char *str)
{
    double x;
    char *p;

    x = strtod(str, &p);

    if      (toupper(*p) == 'K') x *= (1LL << 10);
    else if (toupper(*p) == 'M') x *= (1LL << 20);
    else if (toupper(*p) == 'G') x *= (1LL << 30);

    return static_cast<T>(x + 0.499);
}

#endif
