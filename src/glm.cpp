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

// #include <stdio.h>

#include "glmtlp.hpp"
#include "utils.hpp"
#include "solver.hpp"

// enum class Family
// {
//     Gaussian,
//     Binomial,
//     Poisson
// };

// enum class Method
// {
//     Lasso,
//     RTLP,
//     CTLP
// };

// [[Rcpp::export]]
Rcpp::List call_glm_solver(Rcpp::NumericMatrix &X_mat,
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
                           Rcpp::NumericVector &ncores_val,
                           Rcpp::CharacterVector &family_val,
                           Rcpp::CharacterVector &method_val)
{
    // input
    int family = strcmp(family_val[0], "gaussian") == 0
                     ? 1
                 : strcmp(family_val[0], "binomial") == 0
                     ? 2
                     : 3;
    int method = strcmp(method_val[0], "l1-regularized") == 0
                     ? 1
                 : strcmp(method_val[0], "tlp-regularized") == 0
                     ? 2
                     : 3;

    const int n = X_mat.nrow();
    const int p = X_mat.ncol();
    const int nlambda = lambda_vec.length();
    const int nkappa = kappa_vec.length();

    const double delta = delta_val[0];
    const double tau = tau_val[0];
    const double tol = tol_val[0];

    const int cd_maxit = cd_maxit_val[0];
    const int standardize = standardize_val[0];
    const int ncores = ncores_val[0];

    // const Eigen::Map<Eigen::MatrixXd> X(&X_mat[0], n, p);
    // const Eigen::Map<Eigen::VectorXd> y(&y_vec[0], n);
    // const Eigen::Map<Eigen::VectorXd> w0(&obs_weights[0], n);
    // const Eigen::Map<Eigen::VectorXd> rho0(&pen_factors[0], p);
    // const Eigen::Map<Eigen::VectorXd> lambda(&lambda_vec[0], nlambda);
    // const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], nkappa);

    // output
    // int ntune = (method == Method::CTLP ? nkappa : nlambda);
    int ntune = (method == 3 ? nkappa : nlambda);
    std::vector<Eigen::Triplet<double>> sp_beta_list;
    sp_beta_list.reserve(ntune * std::min(n, p));
    Rcpp::NumericVector intercept(ntune);
    Rcpp::NumericVector deviance(ntune);

    // Eigen::Map<Eigen::VectorXd> b0(&intercept[0], ntune);
    // Eigen::Map<Eigen::VectorXd> dev(&deviance[0], ntune);

    // call solver
    glm_solver(&X_mat[0],
               &y_vec[0],
               &obs_weights[0],
               &pen_factors[0],
               &kappa_vec[0],
               &lambda_vec[0],
               n,
               p,
               nlambda,
               nkappa,
               delta,
               tau,
               tol,
               cd_maxit,
               standardize,
               family,
               method,
               ncores,
               sp_beta_list,
               &intercept[0],
               &deviance[0]);

    Eigen::SparseMatrix<double> beta(p, ntune);
    beta.setFromTriplets(sp_beta_list.begin(), sp_beta_list.end());
    beta.makeCompressed();

    return Rcpp::List::create(
        Rcpp::Named("intercept") = intercept,
        Rcpp::Named("beta") = Rcpp::wrap(beta),
        Rcpp::Named("deviance") = deviance,
        Rcpp::Named("kappa") = kappa_vec,
        Rcpp::Named("lambda") = lambda_vec);
}

// [[Rcpp::export]]
Rcpp::List call_bm_glm_solver(SEXP &X_mat,
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
                              Rcpp::NumericVector &ncores_val,
                              Rcpp::CharacterVector &family_val,
                              Rcpp::CharacterVector &method_val)
{
    // input
    int family = strcmp(family_val[0], "gaussian") == 0
                     ? 1
                 : strcmp(family_val[0], "binomial") == 0
                     ? 2
                     : 3;
    int method = strcmp(method_val[0], "l1-regularized") == 0
                     ? 1
                 : strcmp(method_val[0], "tlp-regularized") == 0
                     ? 2
                     : 3;

    Rcpp::XPtr<BigMatrix> X_(X_mat);
    const int n = X_->nrow();
    const int p = X_->ncol();
    
    const int nlambda = lambda_vec.length();
    const int nkappa = kappa_vec.length();

    const double delta = delta_val[0];
    const double tau = tau_val[0];
    const double tol = tol_val[0];

    const int cd_maxit = cd_maxit_val[0];
    const int standardize = standardize_val[0];
    const int ncores = ncores_val[0];


    // const Eigen::Map<Eigen::MatrixXd> X(&X_mat[0], n, p);
    // const Eigen::Map<Eigen::VectorXd> y(&y_vec[0], n);
    // const Eigen::Map<Eigen::VectorXd> w0(&obs_weights[0], n);
    // const Eigen::Map<Eigen::VectorXd> rho0(&pen_factors[0], p);
    // const Eigen::Map<Eigen::VectorXd> lambda(&lambda_vec[0], nlambda);
    // const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], nkappa);

    // output
    // int ntune = (method == Method::CTLP ? nkappa : nlambda);
    int ntune = (method == 3 ? nkappa : nlambda);
    std::vector<Eigen::Triplet<double>> sp_beta_list;
    sp_beta_list.reserve(ntune * std::min(n, p));
    Rcpp::NumericVector intercept(ntune);
    Rcpp::NumericVector deviance(ntune);

    // Eigen::Map<Eigen::VectorXd> b0(&intercept[0], ntune);
    // Eigen::Map<Eigen::VectorXd> dev(&deviance[0], ntune);

    // call solver
    glm_solver((double *)X_->matrix(),
               &y_vec[0],
               &obs_weights[0],
               &pen_factors[0],
               &kappa_vec[0],
               &lambda_vec[0],
               n,
               p,
               nlambda,
               nkappa,
               delta,
               tau,
               tol,
               cd_maxit,
               standardize,
               family,
               method,
               ncores,
               sp_beta_list,
               &intercept[0],
               &deviance[0]);

    Eigen::SparseMatrix<double> beta(p, ntune);
    beta.setFromTriplets(sp_beta_list.begin(), sp_beta_list.end());
    beta.makeCompressed();

    return Rcpp::List::create(
        Rcpp::Named("intercept") = intercept,
        Rcpp::Named("beta") = Rcpp::wrap(beta),
        Rcpp::Named("deviance") = deviance,
        Rcpp::Named("kappa") = kappa_vec,
        Rcpp::Named("lambda") = lambda_vec);
}