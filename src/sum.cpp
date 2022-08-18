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
#include "solver.h"

// Assume X has standardized columns.
// Assume y is centered.
// Observations are weighted by 1.

// [[Rcpp::export]]
Rcpp::List call_sum_solver(Rcpp::NumericMatrix &XX_mat,
                           Rcpp::NumericVector &Xy_vec,
                           Rcpp::NumericVector &n_val,
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
    const int method = strcmp(method_val[0], "l1-regularized") == 0
                           ? 1
                       : strcmp(method_val[0], "tlp-regularized") == 0
                           ? 2
                           : 3;

    const int n = n_val[0];
    const int p = XX_mat.ncol();
    const int nlambda = lambda_vec.length();
    const int nkappa = kappa_vec.length();

    const double delta = delta_val[0];
    const double tau = tau_val[0];
    const double tol = tol_val[0];
    const int cd_maxit = cd_maxit_val[0];
    // const int ncores = ncores_val[0]; ???

    // const Eigen::Map<Eigen::MatrixXd> XX(&XX_mat[0], p, p);
    // const Eigen::Map<Eigen::VectorXd> Xy(&Xy_vec[0], p);
    // const Eigen::Map<Eigen::VectorXd> rho0(&pen_factors[0], p);
    // const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], p);
    // const Eigen::Map<Eigen::VectorXd> lambda(&lambda_vec[0], nlambda);

    // output
    // int ntune = (method == Method::CTLP ? nkappa : nlambda);
    int ntune = (method == 3 ? nkappa : nlambda);
    std::vector<Eigen::Triplet<double>> sp_beta_list;
    sp_beta_list.reserve(ntune * std::min(n, p));
    Rcpp::NumericVector loss(ntune);

    // Eigen::Map<Eigen::MatrixXd> b(&beta[0], p, ntune); // this should be a column sparse matrix
    // Eigen::Map<Eigen::VectorXd> l(&loss[0], ntune);

    // internal variables
    sum_solver(&XX_mat[0],
               &Xy_vec[0],
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
               method,
               sp_beta_list,
               &loss[0]);

    Eigen::SparseMatrix<double> beta(p, ntune);
    beta.setFromTriplets(sp_beta_list.begin(), sp_beta_list.end());
    beta.makeCompressed();

    return Rcpp::List::create(
        Rcpp::Named("beta") = Rcpp::wrap(beta),
        Rcpp::Named("lambda") = lambda_vec,
        Rcpp::Named("kappa") = kappa_vec,
        Rcpp::Named("loss") = loss);
}

// [[Rcpp::export]]
Rcpp::List call_bm_sum_solver(SEXP &XX_mat,
                           Rcpp::NumericVector &Xy_vec,
                           Rcpp::NumericVector &n_val,
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
    const int method = strcmp(method_val[0], "l1-regularized") == 0
                           ? 1
                       : strcmp(method_val[0], "tlp-regularized") == 0
                           ? 2
                           : 3;

    Rcpp::XPtr<BigMatrix> XX_(XX_mat);
    const int p = XX_->ncol();
    const int n = n_val[0];
    const int nlambda = lambda_vec.length();
    const int nkappa = kappa_vec.length();

    const double delta = delta_val[0];
    const double tau = tau_val[0];
    const double tol = tol_val[0];
    const int cd_maxit = cd_maxit_val[0];
    // const int ncores = ncores_val[0]; ???

    // const Eigen::Map<Eigen::MatrixXd> XX(&XX_mat[0], p, p);
    // const Eigen::Map<Eigen::VectorXd> Xy(&Xy_vec[0], p);
    // const Eigen::Map<Eigen::VectorXd> rho0(&pen_factors[0], p);
    // const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], p);
    // const Eigen::Map<Eigen::VectorXd> lambda(&lambda_vec[0], nlambda);

    // output
    // int ntune = (method == Method::CTLP ? nkappa : nlambda);
    int ntune = (method == 3 ? nkappa : nlambda);
    std::vector<Eigen::Triplet<double>> sp_beta_list;
    sp_beta_list.reserve(ntune * std::min(n, p));
    Rcpp::NumericVector loss(ntune);

    // Eigen::Map<Eigen::MatrixXd> b(&beta[0], p, ntune); // this should be a column sparse matrix
    // Eigen::Map<Eigen::VectorXd> l(&loss[0], ntune);

    // internal variables
    sum_solver((double *)XX_->matrix(),
               &Xy_vec[0],
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
               method,
               sp_beta_list,
               &loss[0]);

    Eigen::SparseMatrix<double> beta(p, ntune);
    beta.setFromTriplets(sp_beta_list.begin(), sp_beta_list.end());
    beta.makeCompressed();

    return Rcpp::List::create(
        Rcpp::Named("beta") = Rcpp::wrap(beta),
        Rcpp::Named("lambda") = lambda_vec,
        Rcpp::Named("kappa") = kappa_vec,
        Rcpp::Named("loss") = loss);
}