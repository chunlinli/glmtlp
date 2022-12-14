# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

call_glm_solver <- function(X_mat, y_vec, obs_weights, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, standardize_val, ncores_val, family_val, method_val) {
    .Call(`_glmtlp_call_glm_solver`, X_mat, y_vec, obs_weights, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, standardize_val, ncores_val, family_val, method_val)
}

call_bm_glm_solver <- function(X_mat, y_vec, obs_weights, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, standardize_val, ncores_val, family_val, method_val) {
    .Call(`_glmtlp_call_bm_glm_solver`, X_mat, y_vec, obs_weights, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, standardize_val, ncores_val, family_val, method_val)
}

call_sum_solver <- function(XX_mat, Xy_vec, n_val, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, method_val) {
    .Call(`_glmtlp_call_sum_solver`, XX_mat, Xy_vec, n_val, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, method_val)
}

call_bm_sum_solver <- function(XX_mat, Xy_vec, n_val, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, method_val) {
    .Call(`_glmtlp_call_bm_sum_solver`, XX_mat, Xy_vec, n_val, pen_factors, kappa_vec, lambda_vec, delta_val, tau_val, tol_val, cd_maxit_val, method_val)
}

