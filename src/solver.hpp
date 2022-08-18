

#pragma once

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "glmtlp.hpp"

void glm_solver(
                double* X_ptr,
                double* y_ptr,
                double* w0_ptr,
                double* rho0_ptr,
                double* kappa_ptr,
                double* lambda_ptr,
                const int n,
                const int p,
                const int nlambda,
                const int nkappa,
                const double delta,
                const double tau,
                const double tol,
                const int cd_maxit,
                const int standardize,
                int family,
                int method,
                int ncores,
                std::vector<Eigen::Triplet<double>> &sp_beta_list,
                double* b0_ptr,
                double* dev_ptr
                );

#endif // SOLVER_HPP
