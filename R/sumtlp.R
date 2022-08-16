# /**********
#     R Interface.
#
#     Copyright (C) 2022 Chunlin Li, Yu Yang, Chong Wu
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program. If not, see <https://www.gnu.org/licenses/>.
# **********/


sumtlp <- function(XX, Xy, nobs,
                   method = c("tlp-constrained", "tlp-regularized", "l1-regularized"),
                   nlambda = ifelse(method == "tlp-constrained", 50, 100),
                   lambda_min_ratio = ifelse(nobs < nvars, 5e-2, 1e-3),
                   lambda = NULL, kappa = NULL,
                   df_max = min(nvars, 40),
                   tau = 0.5 * sqrt(log(nvars) / nobs),
                   delta = 2.0, tol = 1e-4,
                   penalty_factor = rep(1.0, nvars),
                   cd_maxit = 10000, ...) {

    # Coersion
    this_call <- match.call()
    method <- match.arg(method)

    nobs <- as.integer(nobs)
    nvars <- as.integer(length(Xy))

    # check data XX, Xy
    xxdim <- dim(XX)
    if (is.null(xxdim) || any(xxdim != nvars)) {
        stop("X should be a squared matrix")
    }

    if (any(is.na(Xy)) || any(is.na(XX))) {
        stop("Missing data (NA's) detected. Take actions
                (e.g., removing cases, removing features, imputation)
                to eliminate missing data before passing
                XX and Xy to the model")
    }

    # check penalty_factor
    if (method == "tlp-constrained") penalty_factor <- (penalty_factor != 0) * 1
    penalty_factor <- as.double(penalty_factor)
    if (length(penalty_factor) != nvars) {
        stop(paste("the length of penalty_factor (",
            length(penalty_factor), ") not equal to the number of variables (",
            nvars, ")",
            sep = ""
        ))
    }

    ## check/setup lambda and kappa
    if (is.null(lambda)) {
        nlambda <- as.integer(nlambda)
        if (nlambda < 2) stop("nlambda must be at least 2")
        if (lambda_min_ratio >= 1) {
            stop("lambda_min_ratio should be less than 1")
        }
        lambda <- setup_lambda_sum(Xy, lambda_min_ratio, nlambda)
    } else {
        nlambda <- as.integer(length(lambda))
        if (nlambda < 1) stop("the length of input lambda must be at least 1")
        if (!is.double(lambda)) lambda <- as.double(lambda)
        lambda <- sort(lambda, decreasing = TRUE)
        if (lambda[nlambda] < 0) stop("lambdas should be non-negative")
        lambda_max <- max(abs(Xy))
        if (lambda[1] < lambda_max) {
            lambda <- c(lambda_max, lambda)
        } else {
            if (lambda[nlambda] >= lambda_max) {
                message("null model, try smaller lambdas")
                return(get_null_output_sum(
                    this_call,
                    family, method, penalty_factor
                ))
            }
        }
    }
    df_max <- as.integer(df_max)
    if (method == "tlp-constrained") {
        if (is.null(kappa)) {
            if (df_max < 1) stop("max degrees of freedom should be at least 1")
            kappa <- seq(0, df_max)
        } else {
            kappa <- round(kappa)
            kappa <- sort(kappa, decreasing = FALSE)
            if (kappa[1] < 0) stop("kappa should be non-negative")
            if (kappa[1] > 0) kappa <- c(0, kappa)
            nkappa <- length(kappa)
            if (kappa[nkappa] > df_max) {
                message("kappa exceeds max degrees of freedom,
                        df_max is ignored")
            }
        }
    } else {
        kappa <- -1
    }

    ## check tau, delta, tol, and maxiters:
    ## may add more on tau, delta, and tol checking
    tau <- as.double(tau)
    delta <- as.double(delta)
    tol <- as.double(tol)
    cd_maxit <- as.integer(cd_maxit)

    fit <- sum_solver(
        XX, Xy, penalty_factor, kappa, lambda, delta, tau, tol, cd_maxit, method
    )

    varnames <- colnames(XX)
    if (is.null(varnames)) varnames <- paste("V", seq(nvars), sep = "")
    rownames(fit$beta) <- varnames
    if (method == "tlp-constrained") {
        colnames(fit$beta) <- paste(kappa)
    } else {
        colnames(fit$beta) <- lambda_names(lambda)
    }

    out <- structure(
        list(
            beta = fit$beta,
            method = method,
            penalty_factor = penalty_factor,
            loss = fit$loss,
            call = this_call
        ),
        class = "glmtlp"
    )
    if (method == "tlp-constrained") {
        out$kappa <- kappa
        out$lambda <- lambda
        out$tau <- tau
    } else if (method == "tlp-regularized") {
        out$lambda <- lambda
        out$tau <- tau
    } else {
        out$lambda <- lambda
    }
    out
}
