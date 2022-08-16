/**********
    C++ Routines for Linear Algebra Operations.

    Copyright (C) 2020-2021  Chunlin Li

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
**********/

#ifndef _UTILS_H_
#define _UTILS_H_

#include <R.h>
#include <R_ext/BLAS.h>


inline double inner_product_simd(const double *array1,
                                 const double *array2,
                                 const int length,
                                 double init)
{
    int one = 1;
    return init + F77_CALL(ddot)(&length, array1, &one, array2, &one);
}

inline double inner_product_simd(const double *array1,
                                 const double *array2,
                                 const double *array3,
                                 const int length,
                                 double init)
{
    for (int i = 0; i < length; ++i)
    {
        init += array1[i] * array2[i] * array3[i];
    }
    return init;
}

inline void vec_add_simd(const double *array1,
                         const double scalar1,
                         const double scalar2,
                         const int length,
                         double *dest)
{
    
    for (int i = 0; i < length; ++i)
    {
        dest[i] = scalar1 * array1[i] + scalar2;
    }
}

inline void vec_add_simd(const double *array1,
                         const double *array2,
                         const double scalar1,
                         const double scalar2,
                         const int length,
                         double *dest)
{
    
    for (int i = 0; i < length; ++i)
    {
        dest[i] = scalar1 * array1[i] + scalar2 * array2[i];
    }
}



inline double accumulate(const double *array,
                         const int length,
                         double init)
{
    for (int i = 0; i < length; ++i)
    {
        init += array[i];
    }
    return init;
}


// compute soft-thresholding function
inline double soft_thresh(double init, double thresh)
{
    if (init > thresh)
        init -= thresh;
    else if (init < -thresh)
        init += thresh;
    else
        init = 0.0;
    return init;
}

#endif
