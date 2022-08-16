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

inline double link(double mu, Family family)
{
    switch (family)
    {
    case Family::Gaussian:
        return mu;
    case Family::Binomial:
        return std::log(mu / (1.0 - mu));
    case Family::Poisson:
        return std::log(mu); // TODO: check this
    default:
        return NAN;
    }
}

inline double compute_deviance(const Eigen::VectorXd &y,
                               const Eigen::VectorXd &eta,
                               const Eigen::VectorXd &w,
                               Family family)
{
    switch (family)
    {
    case Family::Gaussian:
    {
        return (y - eta).array().square().matrix().dot(w);
    }
    case Family::Binomial:
    {
        Eigen::ArrayXd mu = 1.0 / (1.0 + exp(-eta.array()));
        return -2.0 * w.dot((log(mu) * y.array() + log(1.0 - mu) * (1.0 - y.array())).matrix());
    }
    case Family::Poisson:
    {
        Eigen::ArrayXd mu = exp(eta.array());
        return 2.0 * (w.array() * (log(y.array() / mu + 0.000000001) * y.array())).sum();
    }

    default:
    {
        return NAN;
    }
    }
}

double objective(Eigen::VectorXd &v, double bound, double lambda)
{

    double obj = 0.0;
    for (int j = 0; j < v.size(); ++j)
    {
        double change = std::abs(v(j)) - lambda;
        if (change > 0)
        {
            obj += change;
        }
    }
    return obj - bound;
}

// project_l1

void project_l1(Eigen::VectorXd &v,
                // const Eigen::VectorXi &active_set,
                const std::vector<int> &constrain_vars,
                double bound,
                double tol)
{
    // Eigen::VectorXd v_ = v(constrain_vars);
    Eigen::VectorXd v_(constrain_vars.size());
    for (int k = 0; k < constrain_vars.size(); ++k)
    {
        v_(k) = v(constrain_vars[k]);
    }

    double l1_norm = v_.array().abs().sum();
    if (l1_norm > bound)
    {
        double lambda_upper = v_.array().abs().maxCoeff();
        double lambda_lower = std::max(lambda_upper - bound, 0.0);
        double lambda = (lambda_upper + lambda_lower) / 2.0;

        for (;;)
        {
            if (objective(v_, bound, lambda) < 0)
            {
                lambda_upper = lambda;
            }
            else
            {
                lambda_lower = lambda;
            }

            if (lambda_upper - lambda_lower < tol)
            {
                break;
            }

            lambda = (lambda_upper + lambda_lower) / 2.0;
        }

        for (int k = 0; k < v_.size(); ++k)
        {
            v(constrain_vars[k]) = std::max(std::abs(v_(k)) - lambda, 0.0) * ((v_(k) > 0) - (v_(k) < 0));
        }
    }

    // v(constrain_vars) = v_;
    // return v;
}

void constrain_l1_solver(const Eigen::MatrixXd &X,
                         const Eigen::VectorXd &y,
                         const Eigen::VectorXd &w,
                         const Eigen::VectorXd &rho0, // indicate whether be penalized or not, 0.0 or 1.0
                         const double bound,
                         const double tol,
                         const std::vector<int> &active_set,
                         const double multiplier,
                         const int maxit,
                         double &b0,
                         Eigen::VectorXd &beta)
{
    // input
    const int n = X.rows();
    // const int p = X_mat.ncol();
    const int num_active = active_set.size();
    const int p_ = num_active + 1;

    // const double bound = bound_val[0];
    // const double tol = tol_val[0];
    // const double multiplier = multiplier_val[0];
    // const int maxit = maxit_val[0]; // maximum number of iterations

    // printf("family: %d\n", family);

    // const Eigen::Map<Eigen::MatrixXd> X(&X_mat[0], n, p);
    // const Eigen::Map<Eigen::VectorXd> y(&y_vec[0], n);
    // const Eigen::Map<Eigen::VectorXd> w0(&obs_weights[0], n);
    // const Eigen::Map<Eigen::VectorXi> active_set(&active_set_vec[0], num_active);
    // const Eigen::Map<Eigen::VectorXd> rho0(&pen_factors[0], p);

    Eigen::MatrixXd X_(n, p_);
    for (int k = 0; k < num_active; ++k)
    {
        X_.col(k) = X.col(active_set[k]);
    }
    X_.col(p_ - 1) = Eigen::VectorXd::Ones(n);
    // X_ << X(Eigen::placeholders::all, active_set), Eigen::VectorXd::Ones(n);
    X_ = w.array().sqrt().matrix().asDiagonal() * X_;

    Eigen::VectorXd y_(n);
    y_.array() = w.array().sqrt() * y.array();

    if (bound == 0.0)
    {
        Eigen::VectorXd beta_(p_);
        beta_ = X_.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y_);

        // beta_ =  ((X_.transpose() * X_) / n).solve(X_.transpose() * y_);
        b0 = beta_(p_ - 1);

        for (int k = 0; k < num_active; ++k)
        {
            beta(active_set[k]) = beta_(k);
        }

        return;
    }

    // internal variables
    // Eigen::Map<Eigen::VectorXd> b(&beta[0], p);
    Eigen::VectorXd beta1_(p_);
    for (int k = 0; k < num_active; ++k)
    {
        beta1_(k) = beta(active_set[k]);
    }
    beta1_[p_ - 1] = b0;

    Eigen::VectorXd beta2_(p_);
    beta2_ = beta1_;
    Eigen::VectorXd beta2_old(p_);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(p_);

    std::vector<int> constrain_vars; // index of variables to be constrained
    constrain_vars.reserve(num_active);
    for (int j = 0; j < num_active; ++j)
    {
        if (rho0[active_set[j]])
        {
            constrain_vars.push_back(j);
            // printf("%d\n", (int)constrain_vars.size());
        }
    }

    // Eigen::VectorXd rho_(p_);
    // rho_ << rho0(active_set), 0.0;

    Eigen::MatrixXd XX_inv = ((X_.transpose() * X_) / n + multiplier * Eigen::MatrixXd::Identity(p_, p_)).inverse();

    for (int it = 0; it < maxit; ++it)
    {
        beta1_ = XX_inv * (X_.transpose() * y_ / n + multiplier * (beta2_ - u));
        beta2_old = beta2_;

        /* beginning of beta2_ update */
        beta2_ = beta1_ + u;
        project_l1(beta2_, constrain_vars, bound, tol);
        /* end of beta2_ update */

        u = u + beta1_ - beta2_;

        if ((beta1_ - beta2_).array().abs().maxCoeff() < tol &&
            (beta2_ - beta2_old).array().abs().maxCoeff() < tol)
        {
            // printf("converges");
            break;
        }
    }

    // write output
    b0 = beta2_(p_ - 1);

    for (int k = 0; k < num_active; ++k)
    {
        beta(active_set[k]) = beta2_(k);
    }

    // b(active_set) = beta2_(Eigen::seqN(0, num_active));

    // return Rcpp::List::create(
    //     Rcpp::Named("intercept") = b0,
    //     Rcpp::Named("beta") = beta,
    //     // Rcpp::Named("deviance") = deviance,
    //     Rcpp::Named("bound") = bound_val);
}

// [[Rcpp::export]]
Rcpp::List glm_solver(Rcpp::NumericMatrix &X_mat,
                      Rcpp::NumericVector &y_vec,
                      Rcpp::NumericVector &obs_weights,
                      Rcpp::NumericVector &pen_factors,
                      Rcpp::NumericVector &kappa_vec,
                      Rcpp::NumericVector &lambda_vec,
                      Rcpp::NumericVector &delta_val,
                      Rcpp::NumericVector &tau_val,
                      Rcpp::NumericVector &tol_val,
                      Rcpp::NumericVector &cd_maxit_val,
                      Rcpp::NumericVector &standardize_val,
                      Rcpp::CharacterVector &family_val,
                      Rcpp::CharacterVector &method_val)
{
    // input
    const Family family = strcmp(family_val[0], "gaussian") == 0
                              ? Family::Gaussian
                              : (strcmp(family_val[0], "binomial") == 0
                                     ? Family::Binomial
                                     : Family::Poisson);
    const Method method = strcmp(method_val[0], "l1-regularized") == 0
                              ? Method::Lasso
                          : strcmp(method_val[0], "tlp-regularized") == 0
                              ? Method::RTLP
                              : Method::CTLP;

    const int n = X_mat.nrow();
    const int p = X_mat.ncol();
    const int nlambda = lambda_vec.length();
    const int nkappa = kappa_vec.length();

    const double delta = delta_val[0];
    const double tau = tau_val[0];
    const double tol = tol_val[0];
    // const double multiplier = multiplier_val[0];
    const int cd_maxit = cd_maxit_val[0];
    // const int dc_maxit = dc_maxit_val[0];
    // const int admm_maxit = admm_maxit_val[0];
    const int standardize = standardize_val[0];

    // const int df_max = df_max_val[0];
    // const int user = user_val[0];

    // printf("family: %d\n", family);

    const Eigen::Map<Eigen::MatrixXd> X(&X_mat[0], n, p);
    const Eigen::Map<Eigen::VectorXd> y(&y_vec[0], n);
    const Eigen::Map<Eigen::VectorXd> w0(&obs_weights[0], n);
    const Eigen::Map<Eigen::VectorXd> rho0(&pen_factors[0], p);
    const Eigen::Map<Eigen::VectorXd> lambda(&lambda_vec[0], nlambda);
    const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], nkappa);

    // const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], nkappa);

    // output
    int ntune = (method == Method::CTLP ? nkappa : nlambda);
    Rcpp::NumericMatrix beta(p, ntune);
    Rcpp::NumericVector intercept(ntune);
    Rcpp::NumericVector deviance(ntune);

    Eigen::Map<Eigen::MatrixXd> b(&beta[0], p, ntune);
    Eigen::Map<Eigen::VectorXd> b0(&intercept[0], ntune);
    Eigen::Map<Eigen::VectorXd> dev(&deviance[0], ntune);

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
    std::fill(b.data(), b.data() + p * ntune, 0.0);

    double intercept_new = b0[0];
    double intercept_old = b0[0]; // for warm start
    Eigen::VectorXd beta_new = b.col(0);
    Eigen::VectorXd beta_old = b.col(0); // for warm start

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
        Rcpp::checkUserInterrupt();

        std::copy(rho0.data(), rho0.data() + p, rho.data());
        // warm start (this may be a bit inefficient for Lasso)
        beta_new = beta_old;
        intercept_new = intercept_old;
        rw = rw_old;
        eta = eta_old;

        /* beginning of varaible screening module */
        double cutoff = 2.0 * lambda(k) - lambda(k - 1);
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
                            Rcpp::warning("Maximum number of coordinate descent iterations reached: k = %d.\n", k);
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
                    {
                        EXIT_FLAG++;

                        // TODO: All beta.col(>=k) should be NaN
                        Rcpp::warning("Model saturated: dev(%d) < 0.01 * null_deviance.", k);
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

            for (int kk = 1; kk < nkappa; ++kk)
            {

                int df = kappa(kk);
                if (df > df_max)
                {
                    break;
                }

                // TODO: newton step
                Eigen::VectorXd mu(n);
                Eigen::VectorXd z = y;
                Eigen::VectorXd beta_work = Eigen::VectorXd::Zero(p);

                for (int it_newton = 0; it_newton < 20; ++it_newton)
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
                    }
                    /* end of observation weights update */

                    // TODO: least squares

                    // define matrix X_
                    int p_ = support_size + df + 1;
                    Eigen::MatrixXd X_(n, support_size + df + 1);
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

                    eta.array() = (X_ * beta_work_).array() / (w.array().sqrt()+0.000000001);

                    // update beta
                    for (int jj = 0; jj < support_size + df; ++jj)
                    {
                        beta_work(active_idx[jj]) = beta_work_(jj);
                    }
                    intercept_new = beta_work_(support_size + df);

                    //eta.array() = (X * beta_work).array() + intercept_new;

                    if (family == Family::Gaussian)
                        break;

                    /* check newton-raphson convergence */
                    // if ((b.col(k) - b_old).array().abs().maxCoeff() <= tol)
                    if ((beta_work - b_old).array().abs().maxCoeff() <= tol)
                    {
                        // printf("Newton-Raphson converged.\n");
                        break;
                    }
                }

                double dev_ = compute_deviance(y, eta, w0, family);
                //printf("dev = %f, kappa = %d\n", dev_, df);


                if (dev_ < dev(kk))
                {
                    dev(kk) = dev_;
                    b0(kk) = intercept_new;
                    b.col(kk) = beta_work;
                }
            }

        } /* end of projection module */

        else
        {
            // write result for method == Lasso or RTLP
            dev(k) = compute_deviance(y, eta, w0, family);
            b.col(k) = beta_new;
            b0(k) = intercept_new;
        }

        if (EXIT_FLAG)
        {
            break;
        }

    } /* end of tuning parameter sequence */

    return Rcpp::List::create(
        Rcpp::Named("intercept") = intercept,
        Rcpp::Named("beta") = beta,
        Rcpp::Named("deviance") = deviance,
        Rcpp::Named("kappa") = kappa_vec,
        Rcpp::Named("lambda") = lambda_vec);
}
