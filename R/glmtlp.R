# /**********
#     R Interface: Fit a GLM with Constrained and Regularized TLP/Lasso
#
#     Copyright (C) 2021-2022 Yu Yang, Chunlin Li, Chong Wu
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

#' Fit a GLM with Constrained and Regularized TLP/Lasso
#'
#' @description
#' Fit generalized linear models via penalized maximum likelihood. The
#'   regularization path is computed for the l0, lasso, or truncated lasso
#'   penalty at a grid of values for the regularization parameter \code{lambda}
#'   or \code{kappa}. Fits linear and generalized linear models.
#'
#' @details
#' The sequence of models indexed by \code{lambda} (when \code{penalty = c('l1', 'tlp')})
#'   or \code{kappa} (when \code{penalty = 'l0'}) is fit by the coordinate
#'   descent algorithm.
#'
#'   The objective function for the \code{"gaussian"} family is:
#'   \deqn{1/2 RSS/nobs + \lambda*penalty,} and for the other models it is:
#'   \deqn{-loglik/nobs + \lambda*penalty.}
#'   Also note that, for \code{"gaussian"}, \code{glmtlp} standardizes y to
#'   have unit variance (using 1/(n-1) formula).
#'
#'   ## Details on \code{family} option
#'
#'   \code{glmtlp} currently only supports built-in families, which are specified by a
#'   character string. For all families, the returned object is a regularization
#'   path for fitting the generalized linear regression models, by maximizing the
#'   corresponding penalized log-likelihood. \code{glmtlp(..., family="binomial")}
#'   fits a traditional logistic regression model for the log-odds.
#'
#'   ## Details on \code{penalty} option
#'
#'   The built-in penalties are specified by a character string. For \code{l0}
#'   penalty, \code{kappa} sequence is used for generating the regularization
#'   path, while for \code{l1} and \code{tlp} penalty, \code{lambda} sequence
#'   is used for generating the regularization path.
#'
#' @param X Input matrix, of dimension \code{nobs} x \code{nvars};
#'   each row is  an observation vector.
#' @param y Response variable, of length \code{nobs}. For \code{family="gaussian"},
#'   it should be quantitative; for \code{family="binomial"}, it should be either
#'   a factor with two levels or a binary vector.
#' @param family A character string representing one of the built-in families.
#'   See Details section below.
#' @param penalty A character string representing one of the built-in penalties.
#'   \code{"l0"} represents the \eqn{L_0} penalty, \code{"l1"} represents the
#'   lasso-type penalty (\eqn{L_1} penalty), and \code{"tlp"} represents the
#'   truncated lasso penalty.
#' @param nlambda The number of \code{lambda} values. Default is 100.
#' @param lambda.min.ratio The smallest value for \code{lambda}, as a fraction of
#'   \code{lambda.max}, the smallest value for which all coefficients are zero.
#'   The default depends on the sample size \code{nobs} relative to the number
#'   of variables \code{nvars}. If \code{nobs > nvars}, the default is
#'   \code{0.0001}, and if \code{nobs < nvars}, the default is \code{0.01}.
#' @param lambda A user-supplied \code{lambda} sequence. Typically, users should let
#'   the program compute its own \code{lambda} sequence based on
#'   \code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
#'   \code{lambda} will override this. WARNING: please use this option with care.
#'   \code{glmtlp} relies on warms starts for speed, and it's often faster to
#'   fit a whole path than a single fit. Therefore, provide a decreasing sequence
#'   of \code{lambda} values if you want to use this option. Also, when
#'   \code{penalty = 'l0'}, it is not recommended for the users to supply
#'   this parameter.
#' @param kappa A user-supplied \code{kappa} sequence. Typically, users should
#'   let the program compute its own \code{kappa} sequence based on \code{nvars}
#'   and \code{nobs}. This sequence is used when \code{penalty = 'l0'}.
#' @param tau A tuning parameter used in the TLP-penalized regression models.
#'   Default is  \code{0.3 * sqrt(log(nvars)/nobs)}.
#' @param delta A tuning parameter used in the coordinate majorization descent
#'   algorithm. See Yang, Y., & Zou, H. (2014) in the reference for more detail.
#' @param tol Tolerance level for all iterative optimization algorithms.
#' @param weights Observation weights. Default is 1 for each observation.
#' @param penalty.factor Separate penalty factors applied to each
#'   coefficient, which allows for differential shrinkage. Default is 1
#'   for all variables.
#' @param standardize Logical. Whether or not standardize the input matrix
#'   \code{X}; default is \code{TRUE}.
#' @param dc.maxit Maximum number of iterations for the DC (Difference of
#'   Convex Functions) programming; default is 20.
#' @param cd.maxit Maximum number of iterations for the coordinate descent
#'   algorithm; default is 10^4.
#' @param nr.maxit Maximum number of iterations for the Newton-Raphson method;
#'   default is 500.
#' @param ... Additional arguments.
#' @return An object with S3 class \code{"glmtlp"}.
#'
#' \item{beta}{a \code{nvars x length(kappa)} matrix of
#'   coefficients when \code{penalty = 'l0'}; or a \code{nvars x length(lambda)}
#'   matrix of coefficients when \code{penalty = c('l1', 'tlp')}.}
#' \item{call}{the call that produces this object.}
#' \item{family}{the distribution family used in the model fitting.}
#' \item{intercept}{the intercept vector, of \code{length(kappa)} when
#'   \code{penalty = 'l0'} or \code{length(lambda)} when
#'   \code{penalty = c('l1', 'tlp')}.}
#' \item{lambda}{the actual sequence of \code{lambda} values used. Note that
#'   the length may be smaller than the provided \code{nlambda} due to removal
#'   of saturated values.}
#' \item{penalty}{the penalty type in the model fitting.}
#' \item{penalty.factor}{the penalty factor for each coefficient used in the model fitting.}
#' \item{tau}{the tuning parameter used in the model fitting, available when
#'   \code{penalty = 'tlp'}.}
#'
#' @author Chunlin Li, Yu Yang, Chong Wu
#'   \cr Maintainer: Yu Yang \email{yang6367@umn.edu}
#'
#' @seealso \code{print}, \code{predict}, \code{coef} and \code{plot} methods,
#' and the \code{cv_glmtlp} function.
#'
#' @references Shen, X., Pan, W., & Zhu, Y. (2012).
#'   {Likelihood-based selection and sharp parameter estimation.}
#'   \emph{Journal of the American Statistical Association}, 107(497), 223-232.
#'   \cr Shen, X., Pan, W., Zhu, Y., & Zhou, H. (2013).
#'   {On constrained and regularized high-dimensional regression.}
#'   \emph{Annals of the Institute of Statistical Mathematics}, 65(5), 807-832.
#'   \cr Li, C., Shen, X., & Pan, W. (2021).
#'   {Inference for a large directed graphical model with interventions.}
#'   \emph{arXiv preprint} arXiv:2110.03805.
#'   \cr Yang, Y., & Zou, H. (2014).
#'   {A coordinate majorization descent algorithm for l1 penalized learning.}
#'   \emph{Journal of Statistical Computation and Simulation}, 84(1), 84-95.
#'   \cr Two R packages: \emph{glmnet} and \emph{ncvreg}.
#'
#' @keywords models regression
#'
#' @examples
#'
#' # Gaussian
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' fit1 <- glmtlp(X, y, family = "gaussian", penalty = "l0")
#' fit2 <- glmtlp(X, y, family = "gaussian", penalty = "l1")
#' fit3 <- glmtlp(X, y, family = "gaussian", penalty = "tlp")
#'
#' # Binomial
#'
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- sample(c(0, 1), 100, replace = TRUE)
#' fit <- glmtlp(X, y, family = "binomial", penalty = "l1")
#' @importFrom stats model.matrix
#' @export glmtlp

glmtlp <- function(X, y, family = c("gaussian", "binomial", "poisson"),
                   method = c("tlp-constrained", "tlp-regularized", "l1-regularized"),
                   nlambda = ifelse(method == "tlp-constrained", 50, 100),
                   lambda_min_ratio = ifelse(nobs < nvars, 5e-2, 1e-3),
                   lambda = NULL, kappa = NULL,
                   df_max = min(nvars, 40),
                   tau = 0.5 * sqrt(log(nvars) / nobs),
                   delta = 2.0, tol = 1e-4,
                   weights = NULL, penalty_factor = rep(1.0, nvars),
                   cd_maxit = 10000,
                   standardize = FALSE, eager = TRUE, return_data = FALSE, ...) {
    this_call <- match.call()
    family <- match.arg(family)
    method <- match.arg(method)

    # check data X and y
    xdim <- dim(X)
    if (is.null(xdim) || xdim[2] < 1) {
        stop("X should be a matrix")
    }
    if (!inherits(X, "matrix")) {
        tmp <- try(X <- model.matrix(~ 0 + ., data = X), silent = TRUE)
        if (inherits(tmp, "try-error")) {
            stop("X must be a matrix or able to be coerced to a matrix")
        }
    }
    if (typeof(X) == "character") stop("X must be a numeric matrix")
    if (typeof(X) == "integer") storage.mode(X) <- "double"

    nobs <- as.integer(xdim[1])
    nvars <- as.integer(xdim[2])

    y <- drop(y)
    if (length(y) != nobs) {
        stop(paste("number of observations in y (", length(y),
            ") not equal to the number of rows of X (", nobs, ")",
            sep = ""
        ))
    }
    if (!is.double(y)) {
        op <- options(warn = 2) # turn warnings into errors
        on.exit(options(op))
        y <- tryCatch(
            error = function(cond) {
                stop("y must be numeric or able to be coerced to numeric")
            },
            as.double(y)
        )
        options(op)
    }
    if (any(is.na(y)) || any(is.na(X))) {
        stop("Missing data (NA's) detected. Take actions
                (e.g., removing cases, removing features, imputation)
                to eliminate missing data before passing X and y to the model")
    }
    if (family == "binomial" && length(table(y)) > 2) {
        stop("Attemping to use family='binomial' with non-binary data",
            call. = FALSE
        )
    }
    if (family == "binomial" && !identical(sort(unique(y)), 0:1)) {
        y <- as.double(y == max(y))
    }

    # check penalty_factor and weights
    if (method == "tlp-constrained") penalty_factor <- (penalty_factor != 0) * 1
    penalty_factor <- as.double(penalty_factor)
    if (length(penalty_factor) != nvars) {
        stop(paste("the length of penalty_factor (",
            length(penalty_factor), ") not equal to the number of variables (",
            nvars, ")",
            sep = ""
        ))
    }
    if (is.null(weights)) {
        weights <- rep(1.0, nobs)
    } else if (length(weights) != nobs) {
        stop(paste("number of elements in weights (",
            length(weights), ") not equal to the number of rows of X (",
            nobs, ")",
            sep = ""
        ))
    } else {
        weights <- as.double(weights)
        weights <- weights * nobs / sum(weights)
    }

    ## check/setup lambda and kappa
    if (is.null(lambda)) {
        nlambda <- as.integer(nlambda)
        if (nlambda < 2) stop("nlambda must be at least 2")
        if (lambda_min_ratio >= 1) {
            stop("lambda_min_ratio should be less than 1")
        }
        lambda <- setup_lambda(X, y, weights, lambda_min_ratio, nlambda)
    } else {
        nlambda <- as.integer(length(lambda))
        if (nlambda < 1) stop("the length of input lambda must be at least 1")
        if (!is.double(lambda)) lambda <- as.double(lambda)
        lambda <- sort(lambda, decreasing = TRUE)
        if (lambda[nlambda] < 0) stop("lambdas should be non-negative")
        lambda_max <- get_lambda_max(X, y, weights)
        if (lambda[1] < lambda_max) {
            lambda <- c(lambda_max, lambda)
        } else {
            if (lambda[nlambda] >= lambda_max) {
                message("null model, try smaller lambdas")
                out <- get_null_output(
                    this_call,
                    y, family, method, penalty_factor, weights
                )
                if (return_data) {
                    out$X <- X
                    out$y <- y
                }
                return(out)
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
    standardize <- as.integer(standardize)


    if (eager) {
        fit <- call_glm_solver(
            X, y, weights, penalty_factor, kappa, lambda, delta, tau,
            tol, cd_maxit, standardize, family, method
        )
        ## glm_solver returns intercept, beta, deviance, lambda, kappa

        varnames <- colnames(X)
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
                intercept = fit$intercept,
                family = family,
                method = method,
                penalty_factor = penalty_factor,
                weights = weights,
                deviance = fit$deviance,
                call = this_call,
                is_trained = TRUE
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
    } else {
        is_trained <- FALSE
        out <- get_null_output(
            this_call,
            y, family, method, penalty_factor, weights, is_trained
        )
    }

    if (return_data) {
        out$X <- X
        out$y <- y
    }

    out
}



get_null_output <- function(this_call, y, family, method,
                            penalty_factor, weights, is_trained = TRUE) {
    structure(
        list(
            beta = ifelse(is_trained, 0, NA),
            intercept = ifelse(is_trained, weighted.mean(y, weights), NA),
            # this is wrong for binomial
            family = family,
            method = method,
            penalty_factor = penalty_factor,
            weights = weights,
            call = this_call,
            is_trained = is_trained
        ),
        class = "glmtlp"
    )
}
