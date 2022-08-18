



#include "utils.hpp"

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