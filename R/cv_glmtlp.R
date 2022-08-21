# /**********
#     R Interface: Cross-validation for glmtlp
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

#' Cross-validation for glmtlp
#'
#' Performs k-fold cross-validation for l0, l1, or TLP-penalized regression models
#' over a grid of values for the regularization parameter \code{lambda}
#' (if \code{penalty="l0"}) or \code{kappa} (if \code{penalty="l0"}).
#'
#' The function calls \code{glmtlp} \code{nfolds}+1 times; the first call to get the
#'   \code{lambda} or \code{kappa} sequence, and then the rest to compute
#'   the fit with each of the folds omitted. The cross-validation error is based
#'   on deviance (check here for more details). The error is accumulated over the
#'   folds, and the average error and standard deviation is computed.
#'
#'   When \code{family = "binomial"}, the fold assignment (if not provided by
#'   the user) is generated in a stratified manner, where the ratio of 0/1 outcomes
#'   are the same for each fold.
#'
#' @param X input matrix, of dimension \code{nobs} x \code{nvars}, as in
#'   \code{glmtlp}.
#' @param y response, of length nobs, as in \code{glmtlp}.
#' @param seed the seed for reproduction purposes
#' @param nfolds number of folds; default is 10. The smallest value allowable
#'   is \code{nfolds=3}
#' @param obs_fold an optional vector of values between 1 and \code{nfolds}
#'   identifying what fold each observation is in. If supplied, \code{nfolds} can
#'   be missing.
#' @param ncores number of cores utilized; default is 1. If greater than 1,
#'   then \code{doParallel::foreach} will be used to fit each fold; if equal to
#'   1, then for loop will be used to fit each fold. Users don't have to register
#'   parallel clusters outside.
#' @param \dots Other arguments that can be passed to \code{glmtlp}.
#'
#' @return an object of class \code{"cv.glmtlp"} is returned, which is a list
#'   with the ingredients of the cross-validation fit.
#'
#' \item{call}{the function call}
#' \item{cv.mean}{The mean cross-validated error - a vector of length
#'   \code{length(kappa)} if \code{penalty = "l0"} and \code{length{lambda}}
#'   otherwise.}
#' \item{cv.se}{estimate of standard error of \code{cv.mean}.}
#' \item{fit}{a fitted glmtlp object for the full data.}
#' \item{idx.min}{the index of the \code{lambda} or \code{kappa} sequence that
#'   corresponding to the smallest cv mean error.}
#' \item{kappa}{the values of \code{kappa} used in the fits, available when
#'   \code{penalty = 'l0'}.}
#' \item{kappa.min}{the value of \code{kappa} that gives the minimum
#'   \code{cv.mean}, available when \code{penalty = 'l0'}. }
#' \item{lambda}{the values of \code{lambda} used in the fits.}
#' \item{lambda.min}{value of \code{lambda} that gives minimum \code{cv.mean},
#'   available when penalty is 'l1' or 'tlp'.}
#' \item{null.dev}{null deviance of the model.}
#' \item{obs.fold}{the fold id for each observation used in the CV.}
#'
#' @author Chunlin Li, Yu Yang, Chong Wu
#'   \cr Maintainer: Yu Yang \email{yang6367@umn.edu}
#'
#' @seealso \code{glmtlp} and \code{plot}, \code{predict}, and \code{coef}
#' methods for \code{"cv.glmtlp"} objects.
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
#' cv_fit <- cv_glmtlp(X, y, family = "gaussian", penalty = "l1", seed = 1110)
#'
#' # Binomial
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- sample(c(0, 1), 100, replace = TRUE)
#' cv_fit <- cv_glmtlp(X, y, family = "binomial", penalty = "l1", seed = 1110)
#'
#' @import foreach
#' @importFrom stats sd
#' @export cv_glmtlp

cv_glmtlp <- function(X, y, family = c("gaussian", "binomial", "poisson"),
                      method = c("tlp-constrained", "tlp-regularized", "l1-regularized"),
                      lambda = NULL, kappa = NULL,
                      ..., nfolds = 10, obs_fold = NULL,
                      seed = NULL, ncores = 1) {
    cv_call <- match.call(expand.dots = TRUE)

    fit <- glmtlp(
        X = X, y = y, family = family, method = method,
        lambda = lambda, kappa = kappa,
        ncores = ncores, ...
    ) # ncores = ncores may not be good
    nobs <- nrow(X)
    family <- fit$family
    method <- fit$method
    lambda <- fit$lambda
    kappa <- fit$kappa

    if (!is.null(seed)) set.seed(seed)

    if (family == "binomial" && !identical(sort(unique(y)), 0:1)) {
        y <- as.double(y == max(y))
    }
    if (nfolds > nobs) {
        stop(paste("nfolds (",
            nfolds,
            ") > the number of observations (", nobs, ")",
            sep = ""
        ))
    }

    # generate data fold
    if (is.null(obs_fold)) {
        if (family == "binomial") {
            n0 <- sum(y == 0)
            obs_fold[y == 0] <- sample(rep(1:nfolds, length.out = n0))
            obs_fold[y == 1] <- sample(rep(1:nfolds, length.out = nobs - n0))
        } else {
            obs_fold <- sample(rep(1:nfolds, length.out = nobs))
        }
    } else {
        if (length(obs_fold) != nobs) {
            stop("the length of obs_fold is not equal to nobs")
        }
        obs_fold <- as.integer(factor(obs_fold))
        nfolds <- max(obs_fold)
        if (nfolds == 1) stop("nfolds must be greater than 1")
    }

    # parallel
    if (ncores > 1) {
        doParallel::registerDoParallel(cores = ncores)
        cv_res <- foreach(
            fold = 1:nfolds,
            .combine = "rbind",
            .packages = c("glmtlp")
        ) %dopar% {
            fit_fold <- glmtlp(X, y,
                weights = 1 * (obs_fold != fold),
                lambda = lambda, kappa = kappa,
                family = family, method = method,
                ncores = 1, ...
            )
            yhat <- predict.glmtlp(fit_fold,
                X = X[obs_fold == fold, , drop = FALSE],
                type = "response"
            )
            loss <- loss_glmtlp(y = y[obs_fold == fold], yhat = yhat, family = family)
            loss
        }
    } else {
        cv_res <- c()
        for (fold in 1:nfolds) {
            fit_fold <- glmtlp(X, y,
                weights = 1 * (obs_fold != fold),
                lambda = lambda, kappa = kappa,
                family = family, method = method,
                ncores = 1, ...
            )
            yhat <- predict.glmtlp(fit_fold,
                X = X[obs_fold == fold, , drop = FALSE],
                type = "response"
            )
            loss <- loss_glmtlp(y = y[obs_fold == fold], yhat = yhat, family = family)
            cv_res <- rbind(cv_res, loss)
        }
    }
    cv_mean <- apply(cv_res, 2, mean)
    cv_se <- apply(cv_res, 2, sd)
    idx_min <- which.min(cv_mean)[1]

    out <- structure(list(
        call = cv_call,
        fit = fit,
        obs_fold = obs_fold,
        cv_mean = cv_mean,
        cv_se = cv_se,
        idx_min = idx_min,
        null_dev = loss_glmtlp(y = y, yhat = rep(mean(y), nobs), family = family)
    ),
    class = "cv_glmtlp"
    )
    if (method == "tlp-constrained") {
        out$lambda <- lambda
        out$kappa <- kappa
        out$kappa_min <- kappa[idx_min]
    } else {
        out$lambda <- lambda
        out$lambda_min <- lambda[idx_min]
    }
    out
}
