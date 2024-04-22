#ifndef MY_TIMER_H_
#define MY_TIMER_H_

#include <chrono>

using namespace std;

template <class T>
struct MyTimer
{
    using clock = chrono::high_resolution_clock;
    using duration = chrono::duration<T>;
};

using mytimer = MyTimer<double>;

#endif
