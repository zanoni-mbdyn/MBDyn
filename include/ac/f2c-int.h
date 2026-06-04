/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef AC_F2C_INT_H
#define AC_F2C_INT_H

#ifdef __cplusplus
#warning "<ac/f2c-int.h> is not meant to be included in C++ code"
extern "C" {
#endif /* __cplusplus */

#include <ac/f2c.h>

// #define min(a,b) ((a) < (b) ? (a) : (b))
// #define max(a,b) ((a) > (b) ? (a) : (b))

#include <math.h>

// 1. Safe helper for integers using the bit-twiddling hack
static inline int __min_int(int a, int b) {
    return b + ((a - b) & ((a - b) >> (sizeof(int) * 8 - 1)));
}
static inline long __min_long(long a, long b) {
    return b + ((a - b) & ((a - b) >> (sizeof(long) * 8 - 1)));
}

static inline int __max_int(int a, int b) {
    return a - ((a - b) & ((a - b) >> (sizeof(int) * 8 - 1)));
}
static inline long __max_long(long a, long b) {
    return a - ((a - b) & ((a - b) >> (sizeof(long) * 8 - 1)));
}
// 2. The Master C11 Macro
#define min(a, b) _Generic((a) + (b), \
    float: fminf(a, b),                               \
    double: fmin(a, b),                               \
    long double: fminl(a, b),                         \
    long: __min_long(a, b),                           \
    default: __min_int(a, b)                          \
)

#define max(a, b) _Generic((a) + (b), \
    float: fmaxf(a, b),                               \
    double: fmax(a, b),                               \
    long double: fmaxl(a, b),                         \
    long: __max_long(a, b),                           \
    default: __max_int(a, b)                          \
)


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* AC_F2C_INT_H */

