#pragma once
#include "utils.h"
#include "Eigen/Sparse"

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
    double *dev_ptr);

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
    double *loss_ptr);
