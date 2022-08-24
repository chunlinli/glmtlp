/**********
    C++ Routine for Computing Penalized GLM (dense design version).

    Copyright (C) 2020-2022 Chunlin Li, Yu Yang

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

//#include "glmtlp.hpp"

#include "utils.h"
#include "solver.h"

#include "glmtlp_omp.h"

void glm_solver(
    const double *X_ptr,
    const double *y_ptr,
    const double *w0_ptr,
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
    const int standardize,
    Family family,
    Method method,
    int ncores,
    std::vector<Eigen::Triplet<double>> &sp_beta_list,
    double *b0_ptr,
    double *dev_ptr)
{
    // setup OpenMP
#ifdef GLMTLP_OMP_H_
    if (ncores < 1)
    {
        ncores = omp_get_num_procs();
    }
    omp_set_dynamic(0);
    omp_set_num_threads(ncores);
#endif

    const Eigen::Map<const Eigen::MatrixXd> X(X_ptr, n, p);
    const Eigen::Map<const Eigen::VectorXd> y(y_ptr, n);
    const Eigen::Map<const Eigen::VectorXd> w0(w0_ptr, n);
    const Eigen::Map<const Eigen::VectorXd> rho0(rho0_ptr, p);
    const Eigen::Map<const Eigen::VectorXd> lambda(lambda_ptr, nlambda);
    const Eigen::Map<const Eigen::VectorXd> kappa(kappa_ptr, nkappa);

    // output
    int ntune = (method == Method::CTLP ? nkappa : nlambda);
    // int ntune = (method == 3 ? nkappa : nlambda);
    Eigen::Map<Eigen::VectorXd> b0(b0_ptr, ntune);
    Eigen::Map<Eigen::VectorXd> dev(dev_ptr, ntune);

    // internal variables
    double w_sum = (double)n; // TODO: NOT necessary to be initialized to n in cross-validation

    Eigen::VectorXd b_old(p); // used in newton method
    Eigen::VectorXd rw(n);
    Eigen::VectorXd x_sd(p);
    Eigen::VectorXd xwx(p);
    Eigen::VectorXd rho(p);
    Eigen::VectorXd eta(n); // linear predictor
    Eigen::VectorXd w(n);

    std::copy(w0.data(), w0.data() + n, w.data());

    // initialize
    double weight_sum = w.sum();
    w *= w_sum / weight_sum;
    double y_mean = y.dot(w) / w_sum;

    // TODO: user defined??
    std::fill(b0.data(), b0.data() + ntune, link(y_mean, family));
    // std::fill(b.data(), b.data() + p * ntune, 0.0);

    double intercept_new = b0[0];
    double intercept_old = b0[0]; // for warm start
    // Eigen::VectorXd beta_new = b.col(0);
    // Eigen::VectorXd beta_old = b.col(0); // for warm start
    Eigen::VectorXd beta_new = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd beta_old = Eigen::VectorXd::Zero(p); // for warm start

    rw.array() = (y.array() - y_mean) * w.array();
    double null_deviance = compute_deviance(y, Eigen::VectorXd::Constant(n, b0(0)), w, family);
    std::fill(dev.data(), dev.data() + ntune, null_deviance);

    // TODO: inference??
    std::copy(rho0.data(), rho0.data() + p, rho.data());

    // TODO: need more careful check for standardization
    xwx = X.array().square().matrix().transpose() * w / w_sum;
    if (standardize)
    {
        // TODO
        // should be weighted means and variances
        x_sd = (xwx.array() - (X.transpose() * w / w_sum).array().square()).sqrt();
        rho.array() *= x_sd.array();
    }
    else
    {
        // TODO
        std::fill(x_sd.data(), x_sd.data() + p, 1.0);
    }

    std::fill(eta.data(), eta.data() + n, b0(0));

    // add user initialize
    int it_total_cd = 0;
    // int it_total_dc = 0;
    // int it_total_nr = 0;

    int KKT = 0;
    int CONVERGED = 0;
    int EXIT_FLAG = 0;

    std::vector<int> is_active(p, 0); // active set
    std::vector<int> is_strong(p, 0); // strong set

    Eigen::VectorXd rw_old(n);
    Eigen::VectorXd eta_old(n);

    rw_old = rw;
    eta_old = eta;

    // TODO
    // requires lambda_max to be initialized, k may not start from 1 when user == TRUE
    // if user == TRUE, then add lambda_max in front of user supplied lambda
    for (int k = 1; k < nlambda; ++k)
    {
        check_user_interrupt();

        // std::copy(rho0.data(), rho0.data() + p, rho.data());
        // // // warm start (this may be a bit inefficient for Lasso)
        // beta_new = beta_old;
        // intercept_new = intercept_old;
        // rw = rw_old;
        // eta = eta_old;

        /* beginning of varaible screening module */
        double cutoff = 2.0 * lambda(k) - lambda(k - 1);
#pragma omp parallel for schedule(static)
        for (int j = 0; j < p; ++j)
        {
            if (std::abs(X.col(j).dot(rw) / n) >= cutoff * rho(j))
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

                        /* beginning of observation weights update */
                        switch (family)
                        {
                        case Family::Gaussian:
                            break;

                        case Family::Binomial: // currently only allow w0 = 1
                        {
                            for (int i = 0; i < n; ++i)
                            {
                                if (w(i) > 0.0)
                                {
                                    rw(i) = 1.0 / (1.0 + std::exp(-eta(i)));
                                    if (std::abs(rw(i) - 1.0) < tol)
                                    {
                                        rw(i) = 1.0;
                                        w(i) = tol;
                                    }
                                    else if (std::abs(rw(i)) < tol)
                                    {
                                        rw(i) = 0.0;
                                        w(i) = tol;
                                    }
                                    else
                                    {
                                        w(i) = rw(i) * (1.0 - rw(i));
                                    }
                                    rw(i) = (y(i) - rw(i)) * w0(i); // why w0(i)??
                                }
                            }

                            w_sum = w.sum();
#pragma omp parallel for schedule(static)
                            for (int j = 0; j < p; ++j)
                            {
                                if (is_active[j])
                                {
                                    xwx(j) = (X.col(j).array().square() * w.array()).sum() / n;
                                }
                            }

                            // std::copy(b(k).data(), b.col(k).data() + p, b_old.data());
                            std::copy(beta_new.data(), beta_new.data() + p, b_old.data());
                            break;
                        }

                        case Family::Poisson:
                        {
                            for (int i = 0; i < n; ++i)
                            {
                                if (w(i) > 0.0)
                                {
                                    rw(i) = std::exp(eta(i));
                                    w(i) = rw(i);
                                    rw(i) = (y(i) - rw(i)) * w0(i);
                                }
                            }

                            w_sum = w.sum();
#pragma omp parallel for schedule(static)
                            for (int j = 0; j < p; ++j)
                            {
                                if (is_active[j])
                                {
                                    xwx(j) = (X.col(j).array().square() * w.array()).sum() / n;
                                }
                            }

                            // std::copy(b.col(k).data(), b.col(k).data() + p, b_old.data());
                            std::copy(beta_new.data(), beta_new.data() + p, b_old.data());
                            break;
                        }
                        case Family::Other:
                        {
                            // update w, rw
                            // wls_update();

                            w_sum = w.sum();
#pragma omp parallel for schedule(static)
                            for (int j = 0; j < p; ++j)
                            {
                                if (is_active[j])
                                {
                                    xwx(j) = (X.col(j).array().square() * w.array()).sum() / n;
                                }
                            }

                            // std::copy(b.col(k).data(), b.col(k).data() + p, b_old.data());
                            std::copy(beta_new.data(), beta_new.data() + p, b_old.data());
                            break;
                        }
                        }

                        /* end of observation weights update */

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
                                    // bj_old = b(j, k);
                                    bj_old = beta_new(j);
                                    bj = soft_thresh(
                                        bj_old + X.col(j).dot(rw) / (n * xwx(j) * delta),
                                        lambda(k) * rho(j) / (xwx(j) * delta));
                                    b_change = bj - bj_old;

                                    if (std::abs(b_change) > tol * std::sqrt(null_deviance / n))
                                    // if (b_change)
                                    {
                                        CONVERGED = 0;
                                        // difference = std::abs(b_change) > difference
                                        //                  ? std::abs(b_change)
                                        //                  : difference;

                                        rw.array() -= b_change * w.array() * X.col(j).array();
                                        eta += b_change * X.col(j);

                                        // b(j, k) = bj;
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
                            else
                            {
                                // update intercept
                                b_change = rw.sum() / (w_sum * delta);
                                // b0(k) += b_change;
                                intercept_new += b_change;
                                rw -= w * b_change;
                                eta.array() += b_change;
                            }
                        }

                        if (it_total_cd >= cd_maxit)
                        {
                            EXIT_FLAG++;
                            glmtlp_warning("Maximum number of coordinate descent iterations reached.\n");
                            break;
                        }

                        if (family == Family::Gaussian)
                            break;

                        /* check newton-raphson convergence */
                        // if ((b.col(k) - b_old).array().abs().maxCoeff() <= tol)
                        if ((beta_new - b_old).array().abs().maxCoeff() <= tol)
                        {
                            // printf("Newton-Raphson converged.\n");
                            break;
                        }

                    } /* end of newton-raphson */

                    /* beginning of KKT checking for strong set */
                    KKT = 1;
#pragma omp parallel for schedule(static)
                    for (int j = 0; j < p; ++j)
                    {
                        if (is_strong[j] && !is_active[j] &&
                            std::abs(X.col(j).dot(rw)) / n > lambda(k) * rho(j))
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
#pragma omp parallel for schedule(static)
                for (int j = 0; j < p; ++j)
                {
                    if (!is_strong[j] && !is_active[j] &&
                        std::abs(X.col(j).dot(rw)) / n > lambda(k) * rho(j))
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
                // dev(k) = compute_deviance(y, eta, w0, family);
                if (k != nlambda - 1)
                {

                    /* handling of saturated model */
                    if (family != Family::Gaussian && compute_deviance(y, eta, w0, family) < 0.01 * null_deviance)
                    // if (family != 1 && compute_deviance(y, eta, w0, family) < 0.01 * null_deviance)
                    {
                        EXIT_FLAG++;

                        // TODO: All beta.col(>=k) should be NaN
                        glmtlp_warning("Model saturated: dev < 0.01 * null_deviance.\n");

                        break;
                    }

                    /* beginning of warm start module */
                    // TODO: discuss different methods: Lasso
                    if (it_dc == 0)
                    {
                        // std::copy(b.data() + k * p, b.data() + (k + 1) * p,
                        //           b.data() + (k + 1) * p);
                        // b0(k + 1) = b0(k);

                        beta_old = beta_new;
                        intercept_old = intercept_new;
                        rw_old = rw;
                        eta_old = eta;
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
            // if (method == 1 || EXIT_FLAG)
            {
                break;
            }

            // it_total_dc++;
            it_dc++;

            // update penalty factors
            CONVERGED = 1;
            for (int j = 0; j < p; ++j)
            {
                // CONVERGED = 1;

                if (std::abs(beta_new(j)) * x_sd[j] >= tau)
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
                // printf("%d\n", it_dc);
                break;
            }

            // Update active set ???? Maybe not necessary

        } /* end of difference-of-convex programming */

        /* beginning of projection module */
        if (method == Method::CTLP)
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
                sp_beta_idx[kk] = (int)std::round(sp_beta_idx[kk - 1] + support_size + kappa(kk - 1));
            }
            sp_beta_list.resize((int)std::round(kappa.sum()) + ntune * support_size);

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

            for (int kk = 1; kk < nkappa; ++kk)
            {

                int df = kappa(kk);
                if (df > df_max)
                {
                    df = df_max;
                }

                // TODO: newton step
                Eigen::VectorXd mu = Eigen::VectorXd::Zero(n);
                Eigen::VectorXd z = y;
                Eigen::VectorXd beta_work = Eigen::VectorXd::Zero(p);
                double intercept_work;

                // eta = Eigen::VectorXd::Ones(n) * b0(0);

                for (int it_newton = 0; it_newton < 15; ++it_newton)
                {
                    /* beginning of observation weights update */
                    switch (family)
                    {
                    case Family::Gaussian:
                        break;

                    case Family::Binomial: // currently only allow w0 = 0,1
                    {
                        for (int i = 0; i < n; ++i)
                        {
                            if (w0(i) > 0.0)
                            {
                                mu(i) = 1.0 / (1.0 + std::exp(-eta(i)));
                                if (std::abs(mu(i) - 1.0) < tol)
                                {
                                    mu(i) = 1.0;
                                    w(i) = tol;
                                }
                                else if (std::abs(mu(i)) < tol)
                                {
                                    mu(i) = 0.0;
                                    w(i) = tol;
                                }
                                else
                                {
                                    w(i) = mu(i) * (1.0 - mu(i));
                                }

                                z(i) = eta(i) + (y(i) - mu(i)) / w(i);

                                // rw(i) = (y(i) - rw(i)) * w0(i); // why w0(i)??
                            }
                        }

                        // std::copy(b(k).data(), b.col(k).data() + p, b_old.data());
                        std::copy(beta_work.data(), beta_work.data() + p, b_old.data());
                        break;
                    }

                    case Family::Poisson:
                    {
                        for (int i = 0; i < n; ++i)
                        {
                            if (w(i) > 0.0)
                            {
                                mu(i) = std::exp(eta(i));
                                w(i) = mu(i);
                                z(i) = eta(i) + y(i) / w(i) - 1.0;
                                // rw(i) = (y(i) - rw(i)) * w0(i);
                                // if (std::isnan(w(i)) || std::isinf(w(i))) {
                                //     printf("w(%d) = %f, %d iter\n", i, w(i), it_newton);
                                //     printf("eta(%d) = %f, %d iter\n", i, eta(i), it_newton);
                                //     printf("k = %d", k);
                                // }
                            }
                        }

                        // std::copy(b.col(k).data(), b.col(k).data() + p, b_old.data());
                        std::copy(beta_work.data(), beta_work.data() + p, b_old.data());
                        break;
                    }

                    case Family::Other:
                    {
                        // update

                        // wls_tlp_update();
                        for (int i = 0; i < n; ++i)
                        {
                            if (w(i) > 0.0)
                            {
                                mu(i) = std::exp(eta(i));
                                w(i) = mu(i);
                                z(i) = eta(i) + y(i) / w(i) - 1.0;
                                // rw(i) = (y(i) - rw(i)) * w0(i);
                                // if (std::isnan(w(i)) || std::isinf(w(i))) {
                                //     printf("w(%d) = %f, %d iter\n", i, w(i), it_newton);
                                //     printf("eta(%d) = %f, %d iter\n", i, eta(i), it_newton);
                                //     printf("k = %d", k);
                                // }
                            }
                        }

                        // std::copy(b.col(k).data(), b.col(k).data() + p, b_old.data());
                        std::copy(beta_work.data(), beta_work.data() + p, b_old.data());
                        break;
                    }
                    }
                    /* end of observation weights update */

                    // TODO: least squares

                    // define matrix X_
                    int p_ = support_size + df + 1;
                    Eigen::MatrixXd X_(n, p_);
                    for (int jj = 0; jj < support_size + df; ++jj)
                    {
                        X_.col(jj).array() = X.col(active_idx[jj]).array() * w.array().sqrt();
                    }
                    X_.col(support_size + df) = w.array().sqrt().matrix();

                    Eigen::VectorXd z_(n);
                    z_.array() = z.array() * w.array().sqrt();

                    Eigen::VectorXd beta_work_ = (X_.transpose() * X_ +
                                                  tol * Eigen::MatrixXd::Identity(p_, p_))
                                                     .ldlt()
                                                     .solve(X_.transpose() * z_);

                    eta.array() = (X_ * beta_work_).array() / (w.array().sqrt() + 0.000000001);

                    // update beta
                    for (int jj = 0; jj < support_size + df; ++jj)
                    {
                        beta_work(active_idx[jj]) = beta_work_(jj);
                    }
                    intercept_work = beta_work_(support_size + df);

                    if (family == Family::Gaussian)
                        break;

                    /* check newton-raphson convergence */
                    if ((beta_work - b_old).array().abs().maxCoeff() <= tol)
                    {
                        // printf("Newton-Raphson converged. %d\n", it_newton);
                        break;
                    }
                }

                double dev_ = compute_deviance(y, eta, w0, family);
                // printf("dev = %f, kappa = %d\n", dev_, df);

                if (dev_ + tol < dev(kk))
                {
                    dev(kk) = dev_;
                    b0(kk) = intercept_work;
                    // b.col(kk) = beta_work;
                    int i = 0;
                    for (int j = 0; j < p; ++j)
                    {
                        if (std::abs(beta_work(j)) > tol)
                        {

                            sp_beta_list[sp_beta_idx[kk] + i] = Eigen::Triplet<double>(j, kk, beta_work(j));
                            // sp_beta_list.push_back(Eigen::Triplet<double>(j, kk, beta_work(j)));
                            ++i;
                        }
                    }
                }
            }

        } /* end of projection module */

        else
        {
            // write result for method == Lasso or RTLP
            dev(k) = compute_deviance(y, eta, w0, family);
            b0(k) = intercept_new;
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
