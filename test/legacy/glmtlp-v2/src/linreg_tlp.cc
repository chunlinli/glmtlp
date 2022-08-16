

#include <cmath>    // abs
#include <iostream> // printf
#include "algo.h"
#include "glmtlp.h"
#include "utils.h"

void linreg_tlp_ssr(double *__restrict b0,
                    double *__restrict b,
                    double *__restrict r,
                    const double *__restrict X,
                    const double *__restrict w,
                    const double *__restrict rho,
                    const double *__restrict lambda,
                    const int nlambda,
                    const double tau,
                    const int n,
                    const int p,
                    const double delta,
                    const double tol,
                    const int cd_maxit,
                    const int dc_maxit)
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

    double *r_working = new double[n]; // working residual of d.c.

    // working coordinates
    int *is_working = new int[p];
    int *working_set = new int[p];
    int working_len = 0;

    double *rho_working = new double[p];

    for (int k = 1; k < nlambda; ++k)
    {
        // call lasso solver
        // NEED TO OUTPUT STRONG SET
        linreg_l1_ssr(&b0[k - 1], &b[(k - 1) * p], r, X, v0, v, w, rho,
                      &lambda[k - 1], 2, n, p, delta, tol, cd_maxit);

        // save solutions for warm start
        if (k != nlambda - 1)
        {
            b0[k + 1] = b0[k];
            std::copy(b + k * p, b + k * p + p, b + k * p + p);
        }

        // difference-of-convex program
        std::copy(r, r + n, r_working);
        std::copy(rho, rho + p, rho_working);

        int it = 0;

        for (; it < dc_maxit; ++it)
        {
            int converge = 1; // indicate whether rho changes

            // update rho
            for (int j = 0; j < p; ++j)
            {
                if (std::abs(b[k * p + j]) >= tau)
                {
                    converge = rho_working[j] == 0.0 ? converge : 0;
                    rho_working[j] = 0.0;
                }
                else 
                {
                    converge = rho_working[j] == rho[j] ? converge : 0;
                    rho_working[j] = rho[j];
                }
            }

            // check convergence
            if (converge)
                break;

            // working set
            std::fill(is_working, is_working + p, 0);
            working_len = 0;
            for (int j = 0; j < p; ++j)
            {
                if (b[k * p + j] != 0.0)
                {
                    working_set[working_len++] = j;
                    is_working[j] = 1;
                }
            }

            // iterate on working set
            int it_cd = 0;
            for (;;)
            {
                coordinate_descent(b0 + k, b + k * p, r_working, X, v0, v, w, rho_working,
                                   lambda[k], n, p, delta, tol, cd_maxit,
                                   &it_cd, working_set, working_len);

                // check kkt for all
                int kkt = 1;
                for (int j = 0; j < p; ++j)
                {
                    if (!is_working[j] &&
                        std::abs(inner_product_simd(r_working, X + j * n, n, 0.0)) / n >
                            lambda[k] * rho_working[j])
                    {
                        // KKT condition is violated
                        // add variable to working set
                        kkt = 0;
                        working_set[working_len++] = j;
                        is_working[j] = 1;
                    }
                }

                if (kkt || it_cd >= cd_maxit)
                    break;
            }
        }

        // IF DIVERGES, WARNING?
        
        
        //if (it >= dc_maxit)
        //{
        //printf("it = %d, lambda = %f.\n", it, lambda[k]);
        // printf(
        //     "Warning: the coordinate descent algorithm does not converge (lambda = %f).\n",
        //     lambda[k]);
        //}
    }

    delete[] v;
    delete[] r_working;
    delete[] is_working;
    delete[] working_set;
    delete[] rho_working;
}
