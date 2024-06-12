#ifndef COMMON_H_
#define COMMON_H_

/*
 * Include this file in all executable programs before any variables/functions
 * are declared/defined. The relevant definitions here are the compile-time
 * point dimension parameter POINT_DIM and the point floating-point type Real.
 * These can be adjusted when building the program:
 *
 *    ex1: make DIM=2 FP=32 // points have 2 dimensions and use 32-bit floating point
 *    ex2: make DIM=4 FP=64 // points have 4 dimensions and use 64-bit floating point
 *
 *    etc.
 */

#include <type_traits>
#include "ptraits.h"
#include "timers.h"

#ifndef POINT_DIM
#error "POINT_DIM is undefined"
#else
static_assert(POINT_DIM >= 1, "POINT_DIM must be positive");
#endif

#ifndef FPSIZE
#error "FPSIZE is undefined"
#else
static_assert(FPSIZE == 32 || FPSIZE == 64, "FPSIZE must be 32 or 64");
using Real = std::conditional<(FPSIZE == 32), float, double>::type;
#endif

#endif
