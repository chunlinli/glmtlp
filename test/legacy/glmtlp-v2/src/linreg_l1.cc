/**********
    C++ Routine for Computing Penalized Regression (dense design version).

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

#include <cmath>    // abs
#include <iostream> // printf
#include "algo.h"
#include "glmtlp.h"
#include "utils.h"

/**********
    L1 Gaussian Regression: 
    sum_i [w_i(y_i - b0 - <x_i, b>)^2]/2n + lambda sum_j [rho_j |b_j|]

    Algorithm:
    (1) Set E = S(lambda), where S(lambda) is strong set.
    (2) Solve problem at lambda by coordinate descent 
        using only the features in E.
    (3) Check KKT conditions at this solution for all predictors. 
        (a) If no violation, continue to next lambda with warm start or stop.
        (b) Otherwise add KKT-violated features to E, repeat (2) and (3).

    See <README.md> for details.

    References: 
    <https://www.jstatsoft.org/article/view/v033i01>
    <https://doi.org/10.1111/j.1467-9868.2011.01004.x>
**********/

void linreg_l1_ssr(double *__restrict b0,
                       double *__restrict b,
                       double *__restrict r,
                       const double *__restrict X,
                       const double *__restrict w,
                       const double *__restrict rho,
                       const double *__restrict lambda,
                       const int nlambda,
                       const int n,
                       const int p,
                       const double delta,
                       const double tol,
                       const int cd_maxit)
{
    /*
        b0:         nlambda-dim intercept array;
        b:          (p * nlambda)-dim coefficient array;
        r:          n-dim working residual array;
        X:          (n * p)-dim design matrix array;
        w:          n-dim observation weight array;
        rho:        p-dim penalty factor array;
        lambda:     penalization parameter array;
        nlambda:    number of candidate lambda;
        n:          number of observations;
        p:          number of features;
        delta:      factor for majorized coordinate descent;
        tol:        tolerance error;
        cd_maxit:   maximum iteration of coordinate descent; 
    */
    
    for(int i = 0; i < n; ++i)
    {
        r[i] *= w[i];
    }

    // store v
    // MAY NEED OPTIMIZATION: ONLY VARIABLE IN ACT_SET ARE NEEDED
    double v0 = delta * accumulate(w, n, 0.0);
    double *v = new double[p];
    for (int j = 0; j != p; ++j)
    {
        v[j] = delta * inner_product_simd(
                           X + j * n, X + j * n, w, n, 0.0);
    }

    linreg_l1_ssr(b0, b, r, X, v0, v, w, rho, lambda, nlambda, n, p, delta, tol, cd_maxit);

    delete[] v;
}

void linreg_l1_ssr(double *__restrict b0,
                       double *__restrict b,
                       double *__restrict r,
                       const double *__restrict X,
                       const double v0,
                       const double *__restrict v,
                       const double *__restrict w,
                       const double *__restrict rho,
                       const double *__restrict lambda,
                       const int nlambda,
                       const int n,
                       const int p,
                       const double delta,
                       const double tol,
                       const int cd_maxit)
{
    /*
        b0:         nlambda-dim intercept array;
        b:          (p * nlambda)-dim coefficient array;
        r:          n-dim working residual array;
        X:          (n * p)-dim design matrix array;
        w:          n-dim observation weight array;
        rho:        p-dim penalty factor array;
        lambda:     penalization parameter array;
        nlambda:    number of candidate lambda;
        n:          number of observations;
        p:          number of features;
        delta:      factor for majorized coordinate descent;
        tol:        tolerance error;
        cd_maxit:   maximum iteration of coordinate descent; 
    */

    // active set, strong set, and working set
    int *ever_active = new int[p];
    std::fill(ever_active, ever_active + p, 0);
    int *is_strong = new int[p];
    std::fill(is_strong, is_strong + p, 0);
    int *is_working = new int[p];
    std::fill(is_working, is_working + p, 0);

    // working coordinates
    int *working_set = new int[p];
    int working_len = 0;

    // TODO: NEED TO HANDLE EXCEPTIONS
    // start from k = 1, no need to compute null model (lambda[0])
    for (int k = 1; k < nlambda; ++k)
    {
        
        /* 
            strong rule: 
            |sum_i [x_{ij} r_i(lambda_k)]| < 2lambda_{k+1} - lambda_k
            where r(lambda_k) satisfies optimality condition 
        */
        
        // NEED TO HANDLE EXCEPTIONS (DIVERGENCE CASE)
        for (int j = 0; j < p; ++j)
        {
            // IF RHO > 0.
            if (std::abs(inner_product_simd(r, X + j * n, n, 0.0)) / n >=
                (2.0 * lambda[k] - lambda[k - 1]) * rho[j])
            {
                is_strong[j] = 1;
            }
            else
            {
                is_strong[j] = 0;
            }
        }

        // initialize working set
        working_len = 0;
        if (k > 1)
        {
            for (int j = 0; j < p; ++j)
            {
                if (ever_active[j] || b[k*p + j] != 0.0)
                {
                    ever_active[j] = 1;
                    is_working[j] = 1;
                    working_set[working_len++] = j;
                }
                else
                {
                    is_working[j] = 0;
                }
            }
        }
        else
        {
            for (int j = 0; j < p; ++j)
            {
                if (is_strong[j])
                {
                    working_set[working_len++] = j;
                    is_working[j] = 1;
                }
            }
        }

        int it = 0;  // iteration number
        int kkt = 0; // kkt condition
        for (;;)
        {
            for (;;)
            {

                // lasso solution constrained on strong set
                coordinate_descent(b0 + k, b + k * p, r, X, v0, v, w, rho,
                                   lambda[k], n, p, delta, tol, cd_maxit,
                                   &it, working_set, working_len);

                // check KKT condition for strong set
                kkt = 1;
                for (int j = 0; j < p; ++j)
                {
                    if (is_strong[j] && !is_working[j] &&
                        std::abs(inner_product_simd(r, X + j * n, n, 0.0)) / n >
                            lambda[k] * rho[j])
                    {
                        // KKT condition is violated
                        // add variable to working set
                        kkt = 0;
                        working_set[working_len++] = j;
                        is_working[j] = 1;
                    }
                }

                // termination
                if (kkt || it >= cd_maxit)
                    break;
            }

            // check KKT condition for all
            if (kkt)
            {
                for (int j = 0; j < p; ++j)
                {
                    if (!is_working[j] && !is_strong[j] &&
                        std::abs(inner_product_simd(r, X + j * n, n, 0.0)) / n >
                            lambda[k] * rho[j])
                    {
                        // KKT condition is violated
                        // add variable to working set
                        kkt = 0;
                        working_set[working_len++] = j;
                        is_working[j] = 1;
                    }
                }
            }

            // termination
            if (kkt || it >= cd_maxit)
                break;
        }

        if (kkt)
        {
            if (k != nlambda - 1)
            {
                // estimate converges at lambda[k]
                // warm start
                std::copy(b + k * p, b + k * p + p, b + k * p + p);
                b0[k + 1] = b0[k];
            }
        }
        else
        {
            // NEED TO HANDLE DIVERGENT CASE (START FROM LAST KKT SOLUTION)

            //printf(
            //    "Warning: the coordinate descent algorithm does not converge (lambda = %f).\n",
            //    lambda[k]);
        }
    }

    delete[] ever_active;
    delete[] is_strong;
    delete[] is_working;
    delete[] working_set;
}
