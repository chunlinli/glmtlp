/**********
    C++ Routines for Generalized Linear Model Optimization.

    Copyright (C) 2021  Chunlin Li

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

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <random>
#include "glmtlp.h"
#include "utils.h"

/**********
    Coordinate-wise Descent

    See <README.md> for details.

    References: 
    <https://link.springer.com/article/10.1023/A:1017501703105>
    <https://doi.org/10.1214/07-AOAS131>
**********/


void coordinate_descent(double *__restrict b0,
                        double *__restrict b,
                        double *__restrict r,
                        const double *__restrict X,
                        const double v0,
                        const double *__restrict v,
                        const double *__restrict w,
                        const double *__restrict rho,
                        const double lambda,
                        const int n,
                        const int p,
                        const double delta,
                        const double tol,
                        const int cd_maxit,
                        int *__restrict it,
                        int *__restrict strong,
                        const int strong_len)
{
    double b_j;
    double b_j_working;
    double b_change;

    int iter = *it;
    // coordinate descent iteration
    for (; iter != cd_maxit; ++iter)
    {
        // absolute difference in solution b
        double difference = 0.0;

        // strong set iteration
        // std::shuffle(strong, strong + strong_len,
        //              std::default_random_engine(std::time(0)));
        std::shuffle(strong, strong + strong_len,
                 std::default_random_engine(
                     std::chrono::system_clock::now().time_since_epoch().count()));
        for (int idx = 0; idx != strong_len; ++idx)
        {
            int j = strong[idx];
            b_j = b[j];

            // update b_j
            // NEEDS TO BE REFACTORED!!!
            b_j_working = soft_thresh(
                              inner_product_simd(r, X + j * n, n, b_j * v[j]),
                              n * lambda * rho[j]) /
                          v[j];
            b_change = b_j_working - b_j;

            if (b_change != 0.0)
            {
                difference = std::abs(b_change) > difference ? std::abs(b_change) : difference;


                // update r
                for (int i = 0; i < n; ++i)
                {
                    r[i] -= w[i] * X[j*n + i] * b_change;
                }


                // update r
                //vec_add_simd(r, X + j * n, 1.0, -b_change, n, r);


                b[j] = b_j_working;
            }
        }

        // update intercept
        b_change = accumulate(r, n, 0.0) / v0;
        *b0 += b_change;

        //vec_add_simd(r, 1.0, -b_change, n, r);  //----------------------------------------------
        for (int i = 0; i < n; ++i)
        {
            r[i] -= w[i] * b_change;
        }


        //
        if (difference <= tol)
            break;
    }

    *it = iter;
}


void coordinate_descent(double *__restrict b0,
                        double *__restrict b,
                        double *__restrict r,
                        double *__restrict eta,
                        const double *__restrict X,
                        const double v0,
                        const double *__restrict v,
                        const double *__restrict w,
                        const double *__restrict rho,
                        const double lambda,
                        const int n,
                        const int p,
                        const double delta,
                        const double tol,
                        const int cd_maxit,
                        int *__restrict it,
                        int *__restrict strong,
                        const int strong_len)
{
    double b_j;
    double b_j_working;
    double b_change;

    int iter = *it;
    // coordinate descent iteration
    for (; iter != cd_maxit; ++iter)
    {
        // absolute difference in solution b
        double difference = 0.0;

        // strong set iteration
        // std::shuffle(strong, strong + strong_len,
        //              std::default_random_engine(std::time(0)));
        std::shuffle(strong, strong + strong_len,
                 std::default_random_engine(
                     std::chrono::system_clock::now().time_since_epoch().count()));
        for (int idx = 0; idx != strong_len; ++idx)
        {
            int j = strong[idx];
            b_j = b[j];

            // update b_j
            // NEEDS TO BE REFACTORED!!!
            b_j_working = soft_thresh(
                              inner_product_simd(r, X + j * n, n, b_j * v[j]),
                              n * lambda * rho[j]) /
                          v[j];
            b_change = b_j_working - b_j;

            if (b_change != 0.0)
            {
                difference = std::abs(b_change) > difference ? std::abs(b_change) : difference;




                // update r
                for (int i = 0; i < n; ++i)
                {
                    r[i] -= w[i] * X[j*n + i] * b_change;
                    eta[i] += X[j*n + i] * b_change;
                }




                //vec_add_simd(r, X + j * n, 1.0, -b_change, n, r);  // --------------------------
                //vec_add_simd(eta, X + j * n, 1.0, b_change, n, eta); // ------------------------



                b[j] = b_j_working;
            }
        }

        // update intercept
        b_change = accumulate(r, n, 0.0) / v0;
        *b0 += b_change;



        //vec_add_simd(r, 1.0, -b_change, n, r);  //----------------------------------------------
        for (int i = 0; i < n; ++i)
        {
            r[i] -= w[i] * b_change;
        }




        //
        if (difference <= tol)
            break;
    }

    *it = iter;
}

void newton_raphson(double *__restrict b0,
                    double *__restrict b,
                    double *__restrict r,
                    double *__restrict eta,
                    const double *__restrict y,
                    const double *__restrict X,
                    double *__restrict w,
                    const double *__restrict rho,
                    const double lambda,
                    const int n,
                    const int p,
                    const double delta,
                    const double tol,
                    int *__restrict it,
                    const int nr_maxit,
                    const int cd_maxit,
                    const int *__restrict is_working,
                    int *__restrict working_set,
                    const int working_len)
{
    // newton-raphson iteration
    int iter = *it;
    int cd_it = 0;

    double v0;
    double *v = new double[p];
    double *b_working = new double[p];
    std::copy(b, b + p, b_working);

    for (; iter != nr_maxit; ++iter)
    {

        // update w, r, ...
        for (int i = 0; i < n; ++i)
        {
            r[i] = 1.0 / (1.0 + std::exp(-*b0 - eta[i]));
            if (std::abs(r[i] - 1.0) < tol)
            {
                r[i] = 1.0;
                w[i] = tol;
            }
            else if (std::abs(r[i]) < tol)
            {
                r[i] = 0.0;
                w[i] = tol;
            }
            else
            {
                w[i] = r[i] * (1.0 - r[i]);
            }
            r[i] = y[i] - r[i];
        }




        // update v0, v
        v0 = delta * accumulate(w, n, 0.0);
        for (int j = 0; j < p; ++j)
        {
            if (is_working[j])
            {
                v[j] = delta * inner_product_simd(
                               X + j * n, X + j * n, w, n, 0.0);
            }
        }






        // weighted lasso over a working set
        coordinate_descent(b0, b_working, r, eta, X, v0, v, w, rho,
                           lambda, n, p, delta, tol, cd_maxit,
                           &cd_it, working_set, working_len);

        




        // check convergence
        double difference = 0.0;
        for (int idx = 0; idx < working_len; ++idx)
        {
            int j = working_set[idx];
            double b_change = b_working[j] - b[j];
            difference = std::abs(b_change) > difference
                             ? std::abs(b_change)
                             : difference;
        }

        std::copy(b_working, b_working + p, b);
        if (difference <= tol)
            break;

        
    }

    delete[] v;
    delete[] b_working;
}
