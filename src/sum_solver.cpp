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

// #include <stdio.h>

#include "glmtlp.h"
#include "utils.h"

void sum_solver(
    const double *XX_ptr,
    const double *Xy_ptr,
    const double *rho0_ptr,
    const double *kappa_ptr,
    const double *lambda_ptr,
    const int n,
    const int p,
    const int nlambda,
    const int nkappa,
    const double delta,
    const double tau,
    const double tol,
    const int cd_maxit,
    Method method,
    std::vector<Eigen::Triplet<double>> &sp_beta_list,
    double *loss_ptr)
{

    const Eigen::Map<const Eigen::MatrixXd> XX(XX_ptr, p, p);
    const Eigen::Map<const Eigen::VectorXd> Xy(Xy_ptr, p);
    const Eigen::Map<const Eigen::VectorXd> rho0(rho0_ptr, p);
    const Eigen::Map<const Eigen::VectorXd> lambda(lambda_ptr, nlambda);
    const Eigen::Map<const Eigen::VectorXd> kappa(kappa_ptr, nkappa);

    // output
    int ntune = (method == Method::CTLP ? nkappa : nlambda);
    //int ntune = (method == 3 ? nkappa : nlambda);
    Eigen::Map<Eigen::VectorXd> l(loss_ptr, ntune);
    std::fill(l.data(), l.data() + ntune, 0.0);

    // internal variables

    Eigen::VectorXd Xr = Xy;
    Eigen::VectorXd rho(p);
    std::copy(rho0.data(), rho0.data() + p, rho.data());

    Eigen::VectorXd beta_new = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd beta_old = Eigen::VectorXd::Zero(p); // for warm start

    // add user initialize
    int it_total_cd = 0;

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

        //std::copy(rho0.data(), rho0.data() + p, rho.data());

        //beta_new = beta_old;
        //Xr = Xr_old;

        /* beginning of varaible screening module */
        double cutoff = 2.0 * lambda(k) - lambda(k - 1);
        for (int j = 0; j < p; ++j)
        {
            if (std::abs(Xr(j)) >= cutoff * rho(j))
            {
                is_strong[j] = 1;
            }
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

                        //if (family == 1)// Family::Gaussian)
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
            //if (method == 1 || EXIT_FLAG)
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
        //if (method == 3)
        {
            // TODO

            const int kappa_max = kappa(nkappa - 1);

            std::vector<int> active_idx;
            active_idx.reserve(std::min(p, n));

            for (int j = 0; j < p; ++j)
            {
                if (rho0(j) < tol)
                {
                    active_idx.push_back(j);
                }
            }
            int support_size = active_idx.size();

            std::vector<int> sp_beta_idx(nkappa, 0);
            for (int kk = 1; kk < nkappa; ++kk)
            {
                sp_beta_idx[kk] = (int) std::round(sp_beta_idx[kk-1] + support_size + kappa(kk - 1));
            }
            sp_beta_list.resize((int) std::round(kappa.sum()) + ntune * support_size);

            std::priority_queue<std::pair<double, int>> queue;
            for (int j = 0; j < p; ++j)
            {
                if (std::abs(beta_new(j)) > tol && rho0(j) > tol)
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
                    df = df_max;
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

                    //if (family == Family::Gaussian)
                    break;
                }

                double loss_ = beta_work.dot(XX * beta_work) * 0.5 - Xy.dot(beta_work);
                if (loss_ + tol < l(kk))
                {
                    l(kk) = loss_;
                    // b.col(kk) = beta_work;
                    int i = 0;
                    for (int j = 0; j < p; ++j)
                    {
                        if (std::abs(beta_work(j)) > tol)
                        {
                            sp_beta_list[sp_beta_idx[kk] + i] = Eigen::Triplet<double>(j, kk, beta_work(j));
                            ++i;
                        }
                    }
                }
            }
        }
        else
        {
            // TODO: update loss
            l(k) = beta_new.dot(XX * beta_new) * 0.5 - Xy.dot(beta_new);
            // b.col(k) = beta_new;
            for (int j = 0; j < p; ++j)
            {
                if (std::abs(beta_new(j)) > tol)
                {
                    sp_beta_list.push_back(Eigen::Triplet<double>(j, k, beta_new(j)));
                }
            }
        }

        if (EXIT_FLAG)
        {
            break;
        }

    } /* end of tuning parameter sequence */

}
