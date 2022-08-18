/**********
    C++ Routine for Estimation of Gaussian Graph Model.

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

#include "glmtlp.hpp"
#include "utils.hpp"



void soft_thresh(Eigen::MatrixXd &A, Eigen::MatrixXd &Thr)
{
    for (int i = 0; i < A.rows(); ++i)
    {
        for (int j = 0; j < A.cols(); ++j)
        {
            if (i != j)
            {
                if (A(i, j) > Thr(i, j))
                    A(i, j) -= Thr(i, j);
                else if (A(i, j) < -Thr(i, j))
                    A(i, j) += Thr(i, j);
                else
                    A(i, j) = 0;
            }
        }
    }
}




// [[Rcpp::export]]
Rcpp::List ggm_solver(Rcpp::NumericMatrix &XX_mat,
                      Rcpp::NumericMatrix &rho_mat, // penalty factors
                      Rcpp::NumericVector &lambda_vec,
                      Rcpp::NumericVector &kappa_vec,
                      Rcpp::NumericVector &multiplier_val,
                      Rcpp::NumericVector &maxit_val,
                      Rcpp::NumericVector &tol_val,
                      Rcpp::CharacterVector &method_val)
{


    const double multiplier = multiplier_val[0];
    const double tol = tol_val[0];
    const int maxit = maxit_val[0];
    const int nlambda = lambda_vec.length();
    const int nkappa = kappa_vec.length();

    const Eigen::Map<Eigen::MatrixXd> XX(&XX_mat[0], XX_mat.nrow(), XX_mat.ncol());
    const Eigen::Map<Eigen::MatrixXd> rho0(&rho_mat[0], rho_mat.nrow(), rho_mat.ncol());
    const Eigen::Map<Eigen::VectorXd> lambda(&lambda_vec[0], nlambda);
    const Eigen::Map<Eigen::VectorXd> kappa(&kappa_vec[0], nkappa);



    const Method method = strcmp(method_val[0], "unconstrained-l1") == 0
                              ? Method::Lasso
                          : strcmp(method_val[0], "unconstrained-tlp") == 0
                              ? Method::RTLP
                              : Method::CTLP;


    const int ntune = method == Method::CTLP ? nkappa : nlambda;
    Rcpp::List Omegas(ntune);


    // need more careful initialization
    Eigen::MatrixXd Omega1 = (XX + tol * Eigen::MatrixXd::Identity(XX.rows(), XX.cols())).inverse();
    Eigen::MatrixXd Omega2 = Omega1;
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(XX.rows(), XX.rows());

    Eigen::MatrixXd rho(XX.rows(), XX.rows());


    int it = 0;
    int EXIT_FLAG = 0;

    Eigen::MatrixXd Omega_old(XX.rows(), XX.rows());

    for (int k = 1; k < nlambda; ++k)
    {

        rho.array() = rho0.array() * lambda(k) / multiplier;

        // admm loop begins
        for (;;)
        {
            // update Omega1
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(multiplier * (Omega2 - U) - XX);
            Eigen::VectorXd eigenvalues = es.eigenvalues();
            eigenvalues = ( eigenvalues.array() + (eigenvalues.array().square() + 4 * multiplier).sqrt() )/ (2 * multiplier);
            Omega1 = es.eigenvectors() * eigenvalues.asDiagonal() * es.eigenvectors().transpose();

            // update Omega2
            Omega_old = Omega2;
            Omega2 = Omega1 + U;
            soft_thresh(Omega2, rho);

            // update U
            U = U + Omega1 - Omega2;

            // check convergence
            if ((Omega1 - Omega2).array().abs().maxCoeff() < tol &&
                (Omega2 - Omega_old).array().abs().maxCoeff() < tol)
            {
                break;
            }

            // check convergence
            if (it > maxit)
            {
                ++EXIT_FLAG;
                break;
            }


        }
        // admm loop ends


        // TODO: should store in sparse matrix
        Omegas[k] = Rcpp::wrap(Omega2);

        if (EXIT_FLAG)
        {
            break;
        }


    }




    return Rcpp::List::create(Rcpp::Named("Omega") = Omegas,
                              Rcpp::Named("lambda") = lambda,
                              Rcpp::Named("kappa") = kappa);

}







