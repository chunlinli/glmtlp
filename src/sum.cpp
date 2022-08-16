/**********
    C++ Routine for Regression with Summary Data (dense design version).

    Copyright (C) 2022 Chunlin Li

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

#include <stdio.h>

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <string.h>
#include <cmath>
#include <queue>
// [[Rcpp::depends(RcppEigen)]]

enum class Family
{
    Gaussian,
    Binomial,
    Poisson
};

enum class Method
{
    Lasso,
    RTLP,
    CTLP
};

// optimizer

// template

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

// inline double link(double mu, Family family)
// {
//     switch (family)
//     {
//     case Family::Gaussian:
//         return mu;
//     case Family::Binomial:
//         return NAN;
//     case Family::Poisson:
//         return NAN;
//     default:
//         return NAN;
//     }
// }

// inline double compute_deviance(const Eigen::VectorXd &y,
//                                const Eigen::VectorXd &eta,
//                                const Eigen::VectorXd &w,
//                                Family family)
// {
//     switch (family)
//     {
//     case Family::Gaussian:
//     {
//         return (y - eta).array().square().matrix().dot(w);
//     }
//     case Family::Binomial:
//     {
//         return NAN;
//     }
//     case Family::Poisson:
//     {
//         return NAN;
//     }

//     default:
//     {
//         return NAN;
//     }
//     }
// }

// Assume X has standardized columns.
// Assume y is centered.
// Observations are weighted by 1.

// [[Rcpp::export]]
Rcpp::List sum_solver(Rcpp::NumericMatrix &XX_mat,
                      Rcpp::NumericVector &Xy_vec,
                      Rcpp::NumericVector &pen_factors,
                      Rcpp::NumericVector &kappa_vec,
                      Rcpp::NumericVector &lambda_vec,
                      Rcpp::NumericVector &delta_val,
                      Rcpp::NumericVector &tau_val,
                      Rcpp::NumericVector &tol_val,
                      Rcpp::NumericVector &cd_maxit_val,
                      Rcpp::CharacterVector &method_val)
{
    // input
    // const int n = X_mat.nrow();
    // 
    const int p = XX_mat.ncol();
    const int nlambda = lambda_vec.length();
    const int nkappa = kappa_vec.length();
    const double delta = delta_val[0];
    const double tau = tau_val[0];
    const double tol = tol_val[0];
    const int cd_maxit = cd_maxit_val[0];
    // const int dc_maxit = dc_maxit_val[0];
    // const int standardize = standardize_val[0];

    // const int df_max = df_max_val[0];
    // const int user = user_val[0];

    const Family family = Family::Gaussian;
    const Method method = strcmp(method_val[0], "l1-regularized") == 0
                              ? Method::Lasso
                          : strcmp(method_val[0], "tlp-regularized") == 0
                              ? Method::RTLP
                              : Method::CTLP;

    // printf("family: %d\n", family);

    const Eigen::Map<Eigen::MatrixXd> XX(&XX_mat[0], p, p);
    const Eigen::Map<Eigen::VectorXd> Xy(&Xy_vec[0], p);
    const Eigen::Map<Eigen::VectorXd> rho0(&pen_factors[0], p);
    const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], p);
    const Eigen::Map<Eigen::VectorXd> lambda(&lambda_vec[0], nlambda);

    // output
    int ntune = (method == Method::CTLP ? nkappa : nlambda);
    Rcpp::NumericMatrix beta(p, ntune);
    Rcpp::NumericVector loss(ntune);

    Eigen::Map<Eigen::MatrixXd> b(&beta[0], p, ntune); // this should be a column sparse matrix
    Eigen::Map<Eigen::VectorXd> l(&loss[0], ntune);

    // internal variables

    // Eigen::VectorXd b_old(p);
    // Eigen::VectorXd r(n);
    // Eigen::VectorXd x_sd(p);
    // Eigen::VectorXd xx(p);
    Eigen::VectorXd Xr = Xy;
    Eigen::VectorXd rho(p);
    // Eigen::VectorXd eta(n); // linear predictor

    // initialize
    // double weight_sum = w.sum();
    // w *= w_sum / weight_sum;
    // double y_mean = y.dot(w) / w_sum;

    // TODO: user defined??
    // std::fill(b0.data(), b0.data() + nlambda, link(y_mean, family));

    // need change
    std::fill(b.data(), b.data() + p * ntune, 0.0);

    Eigen::VectorXd beta_new = b.col(0);
    Eigen::VectorXd beta_old = b.col(0); // for warm start

    // rw.array() = (y.array() - y_mean) * w.array();
    // double null_deviance
    std::fill(l.data(), l.data() + ntune, 0.0);

    // TODO: inference??
    std::copy(rho0.data(), rho0.data() + p, rho.data());

    // TODO: need more careful check for standardization
    // xwx = X.array().square().matrix().transpose() * w / w_sum;
    // if (standardize)
    // {
    //     // TODO
    //     // should be weighted means and variances
    //     x_sd = (xwx.array() - (X.transpose() * w / w_sum).array().square()).sqrt();
    //     rho.array() *= x_sd.array();
    // }
    // else
    // {
    //     // TODO
    //     std::fill(x_sd.data(), x_sd.data() + p, 1.0);
    // }

    // std::fill(eta.data(), eta.data() + n, b0(0));

    // add user initialize
    int it_total_cd = 0;
    // int it_total_dc = 0;
    // int it_total_nr = 0;

    int KKT = 0;
    int CONVERGED = 0;
    int EXIT_FLAG = 0;

    std::vector<int> is_active(p, 0); // active set
    std::vector<int> is_strong(p, 0); // strong set

    Eigen::VectorXd Xr_old(p); 
    Xr_old = Xr;

    // TODO
    // requires lambda_max to be initialized, k may not start from 1 when user == TRUE
    // if user == TRUE, then add lambda_max in front of user supplied lambda
    for (int k = 1; k < nlambda; ++k)
    {
        Rcpp::checkUserInterrupt();

        std::copy(rho0.data(), rho0.data() + p, rho.data());

        beta_new = beta_old;
        Xr = Xr_old;

        /* beginning of varaible screening module */
        double cutoff = 2.0 * lambda(k) - lambda(k - 1);
        for (int j = 0; j < p; ++j)
        {
            if (std::abs(Xr(j)) >= cutoff * rho(j))
            {
                is_strong[j] = 1;
            }
            // else
            // {
            //     is_strong[j] = 0;
            // }
        }

        if (k == 1)
        {
            std::copy(is_strong.begin(), is_strong.end(), is_active.begin());
        }
        /* end of variable screening module */

        int it_dc = 0;
        for (;;) /* beginning of difference-of-convex programming */
        {

            for (;;)
            {
                for (;;)
                {
                    for (;;) /* beginning of newton-raphson */
                    {

                        for (; it_total_cd != cd_maxit;)
                        {

                            double bj;
                            double bj_old;
                            double b_change;

                            // double difference = 0.0;
                            CONVERGED = 1;

                            // update coordinates
                            for (int j = 0; j < p; ++j)
                            {
                                if (is_active[j])
                                {

                                    // TODO: replace this
                                    bj_old = beta_new(j);
                                    bj = soft_thresh(
                                        bj_old + Xr(j) / (XX(j, j) * delta),
                                        lambda(k) * rho(j) / (XX(j, j) * delta));
                                    b_change = bj - bj_old;

                                    if (std::abs(b_change) > tol)
                                    {
                                        CONVERGED = 0;

                                        Xr -= XX.col(j) * b_change;

                                        // rw.array() -= b_change * w.array() * X.col(j).array();
                                        // eta += b_change * X.col(j);
                                        beta_new(j) = bj;
                                    }
                                }
                            }

                            it_total_cd++;

                            /* check coordinate descent convergence */
                            if (CONVERGED)
                            {
                                // printf("(%d, %d)\n", it_total_cd, k);
                                break;
                            }
                        }

                        if (it_total_cd >= cd_maxit)
                        {
                            EXIT_FLAG++;
                            Rcpp::warning("Maximum number of coordinate descent iterations reached: k = %d.\n", k);
                            break;
                        }

                        if (family == Family::Gaussian)
                            break;

                        /* check newton-raphson convergence */
                        // if ((b.col(k) - b_old).array().abs().maxCoeff() <= tol)
                        // {
                        //     // printf("Newton-Raphson converged.\n");
                        //     break;
                        // }

                    } /* end of newton-raphson */

                    /* beginning of KKT checking for strong set */
                    KKT = 1;
                    for (int j = 0; j < p; ++j)
                    {

                        // TODO: replace this

                        if (is_strong[j] && !is_active[j] &&
                            std::abs(Xr(j)) > lambda(k) * rho(j))
                        {
                            KKT = 0;
                            is_active[j] = 1;
                        }
                    }

                    if (KKT || EXIT_FLAG)
                    {
                        break;
                    }
                    /* end of KKT checking for strong set */
                }

                /* beginning of KKT checking for the rest */
                KKT = 1;
                for (int j = 0; j < p; ++j)
                {

                    // TODO: replace this

                    if (!is_strong[j] && !is_active[j] &&
                        std::abs(Xr(j)) > lambda(k) * rho(j))
                    {
                        KKT = 0;
                        is_active[j] = 1;
                    }
                }

                if (KKT || EXIT_FLAG)
                {
                    break;
                }
                /* end of KKT checking for the rest */
            }

            if (KKT)
            {
                // printf("KKT condition satisfied\n");

                if (k != nlambda - 1)
                {
                    /* handling of saturated model */
                    // if (family != Family::Gaussian && dev(k) < 0.01 * null_deviance)
                    // {
                    //     EXIT_FLAG++;

                    //     // TODO: All beta.col(>=k) should be NaN
                    //     Rcpp::warning("Model saturated: dev(%d) < 0.01 * null_deviance.", k);
                    //     break;
                    // }

                    /* beginning of warm start module */
                    // TODO: discuss different methods: Lasso, TLP
                    if (it_dc == 0)
                    {
                        beta_old = beta_new;
                        Xr_old = Xr;

                        // std::copy(b.data() + k * p, b.data() + (k + 1) * p,
                        //           b.data() + (k + 1) * p);
                        // b0(k + 1) = b0(k);
                    }
                    /* end of warm start module */
                }
            }
            else
            {
                // TODO: All beta.col(>=k) should be NaN

                break;
            }

            if (method == Method::Lasso || EXIT_FLAG)
            {
                break;
            }

            // it_total_dc++;
            it_dc++;

            // update penalty factors
            CONVERGED = 1;
            for (int j = 0; j < p; ++j)
            {
                if (std::abs(beta_new(j)) >= tau)
                {
                    CONVERGED = rho(j) == 0.0 ? CONVERGED : 0;
                    rho(j) = 0.0;
                }
                else
                {
                    CONVERGED = rho(j) == rho0(j) ? CONVERGED : 0;
                    rho(j) = rho0(j);
                }
            }

            // check difference of convex programming convergence
            if (CONVERGED)
            {
                break;
            }

            // Update active set ???? Maybe not necessary

        } /* end of difference-of-convex programming */

        if (method == Method::CTLP)
        {
            const int kappa_max = kappa(nkappa - 1);

            std::vector<int> active_idx;
            active_idx.reserve(p);

            for (int j = 0; j < p; ++j)
            {
                if (rho0(j) == 0.0)
                {
                    active_idx.push_back(j);
                }
            }
            int support_size = active_idx.size();

            std::priority_queue<std::pair<double, int>> queue;
            for (int j = 0; j < p; ++j)
            {
                if (std::abs(beta_new(j)) > tol && rho0(j) != 0.0)
                {
                    queue.push(std::pair<double, int>(std::abs(beta_new(j)), j));
                }
            }

            int df_max = kappa_max < (int)queue.size() ? kappa_max : (int)queue.size();

            for (int df = 0; df < df_max; ++df)
            {
                active_idx.push_back(queue.top().second);
                queue.pop();
            }

            for (int kk = 0; kk < nkappa; ++kk)
            {

                int df = kappa(kk);
                if (df > df_max)
                {
                    break;
                }

                // TODO: newton step
                //Eigen::VectorXd mu(n);
                Eigen::VectorXd beta_work = Eigen::VectorXd::Zero(p);

                for (;;)
                {
                    // define matrix XX_
                    int p_ = support_size + df;
                    Eigen::MatrixXd XX_(p_, p_);

                    for (int jj1 = 0; jj1 < p_; ++jj1)
                    {
                        for (int jj2 = 0; jj2 < p_; ++jj2)
                        {
                            XX_(jj2, jj1) = XX(active_idx[jj2], active_idx[jj1]);
                        }
                        //X_.col(jj).array() = X.col(active_idx[jj]).array() * w.array().sqrt();
                    }
                    //X_.col(support_size + df) = Eigen::VectorXd::Ones(n);

                    Eigen::VectorXd Xy_(p_);
                    for (int jj = 0; jj < p_; ++jj)
                    {
                        Xy_(jj) = Xy(active_idx[jj]);
                    }
                    // y_.array() = y.array() * w.array().sqrt();
                    Eigen::VectorXd beta_work_ = XX_.ldlt().solve(Xy_);

                    //eta = X_ * beta_work_;

                    // update beta
                    for (int jj = 0; jj < p_; ++jj)
                    {
                        beta_work(active_idx[jj]) = beta_work_(jj);
                    }
                    //intercept_new = beta_work_(support_size + df);

                    if (family == Family::Gaussian)
                        break;

                    /* check newton-raphson convergence */
                    // if ((b.col(k) - b_old).array().abs().maxCoeff() <= tol)
                    // if ((beta_work - b_old).array().abs().maxCoeff() <= tol)
                    // {
                    //     // printf("Newton-Raphson converged.\n");
                    //     break;
                    // }
                }

                double loss_ = beta_work.dot(XX * beta_work) * 0.5 - Xy.dot(beta_work);
                if (loss_ < l(kk))
                {
                    l(kk) = loss_;
                    b.col(kk) = beta_work;
                }
            }
        }
        else
        {
            // TODO: update loss
            l(k) = beta_new.dot(XX * beta_new) * 0.5 - Xy.dot(beta_new);
            b.col(k) = beta_new;
        }

        if (EXIT_FLAG)
        {
            break;
        }

    } /* end of tuning parameter sequence */

    return Rcpp::List::create(
        Rcpp::Named("beta") = beta,
        Rcpp::Named("lambda") = lambda_vec,
        Rcpp::Named("kappa") = kappa,
        Rcpp::Named("loss") = loss);
}
