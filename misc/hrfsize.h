#ifndef HRFSIZE_H_
#define HRFSIZE_H_

#include <filesystem>
#include <sstream>
#include <string>

using namespace std;

struct HumanReadable
{
    uintmax_t size;

    HumanReadable(const char *fname) : size(filesystem::file_size(filesystem::path(fname))) {}

    static string str(const char *fname)
    {
        stringstream ss;
        ss << HumanReadable(fname);
        return ss.str();
    }

    template <class c, class t>
    friend basic_ostream<c,t>& operator<<(basic_ostream<c,t>& os, HumanReadable hr)
    {
        int i{};
        double mantissa = hr.size;
        for (; mantissa >= 1024.0; mantissa /= 1024.0, ++i);
        os << ceil(mantissa * 10.) / 10. << i["BKMGTPE"];
        return i? os << "B" : os;
    }
};



#endif
