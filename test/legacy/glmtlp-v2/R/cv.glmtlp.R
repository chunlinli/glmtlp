#' Cross-validation for glmtlp
#'
#' Does k-fold cross-validation for glmtlp, produces a plot, and returns a
#' value for \code{lambda} (or \code{kappa} if \code{penalty="l0"}).
#'
#' The function runs \code{glmtlp} \code{nfolds}+1 times; the first to get the
#'   \code{lambda} or \code{kappa} sequence, and then the remainder to compute 
#'   the fit with each of the folds omitted. The error is accumulated, and the 
#'   average (haven't finished yet) error and standard deviation over the folds 
#'   is computed. 
#'
#' @param X input matrix, of dimension \code{nobs} x \code{nvars}, as in 
#'   \code{glmtlp}.
#' @param y response, of length nobs, as in \code{glmtlp}.
#' @param seed the seed for reproduction purposes
#' @param nfolds number of folds; default is 10. Although \code{nfolds} can be
#'   as large as the sample size (leave-one-out CV), it is not recommended for
#'   large datasets. The smallest value allowable is \code{nfolds=3}
#' @param obs.fold an optional vector of values between 1 and \code{nfolds}
#'   identifying what fold each observation is in. If supplied, \code{nfolds} can
#'   be missing.
#' @param ncores number of cores utilized; default is 1. If greater than 1, 
#'   then use \code{doParallel::foreach} to fit each fold; if equal to 1, then 
#'   use for loop to fit each fold. No need to register parallel outside.
#' @param \dots Other arguments that can be passed to \code{glmtlp}
#' 
#' @return an object of class \code{"cv.glmtlp"} is returned, which is a list
#'   with the ingredients of the cross-validation fit. 
#' 
#' \item{call}{the function call}
#' \item{cv.mean}{The mean cross-validated error - a vector of length
#'   \code{length(kappa)} if \code{penalty = "l0"} and a vector of length 
#'   \code{length{lambda}} otherwise.} 
#' \item{cv.se}{estimate of standard error of \code{cv.mean}.} 
#' \item{fit}{a fitted glmtlp object for the full data.} 
#' \item{idx.min}{the index of the \code{lambda} or \code{kappa} sequence that 
#'   corresponding to the smallest cv mean error.}
#' \item{kappa}{the values of \code{kappa} used in the fits, available when 
#'   penalty = 'l0'.}
#' \item{kappa.min}{the value of \code{kappa} that gives the minimum 
#'   \code{cv.mean}, available when penalty is 'l0'. }
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
#'   \emph{Likelihood-based selection and sharp parameter estimation. 
#'   Journal of the American Statistical Association, 107(497), 223-232.}
#'   \cr Yang, Y., & Zou, H. (2014). \emph{A coordinate majorization descent algorithm 
#'   for l1 penalized learning. Journal of Statistical Computation and 
#'   Simulation, 84(1), 84-95.}
#'   \cr Two R package Github: \emph{ncvreg} and \emph{glmnet}.
#' 
#' @keywords models regression
#' 
#' @examples
#' 
#' # Gaussian
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' cv.fit <- cv.glmtlp(X, y, family = "gaussian", penalty = "l1", seed=2021)
#'
#' # Binomial
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- sample(c(0,1), 100, replace = TRUE)
#' cv.fit <- cv.glmtlp(X, y, family = "binomial", penalty = "l1", seed=2021)
#' 
#' @import foreach
#' @import parallel
#' @importFrom stats sd
#' @export cv.glmtlp

cv.glmtlp <- function(X, y, ..., seed=NULL, nfolds=10, obs.fold=NULL, ncores=1) {
    cv.call <- match.call(expand.dots = TRUE)
    
    fit <- glmtlp(X = X, y = y, ...)  # will do the error checking and generate lambda, kappa sequences
    nobs <- nrow(X)
    family <- fit$family
    penalty <- fit$penalty
    lambda <- fit$lambda
    kappa <- fit$kappa

    if (family == "binomial" & !identical(sort(unique(y)), 0:1)) {
        y <- as.double(y == max(y))
    }
    if (!is.null(seed)) set.seed(seed)
    if (nfolds > nobs) stop(paste("nfolds (", nfolds, ") cannot be larger than the number of observations (", nobs, ")", sep = ""))

    # generate fold
    if (is.null(obs.fold)) {
        if (family == "binomial") {  # stratified sampling for binary data
            n0 <- sum(y == 0)
            obs.fold0 <- 1:n0 %% nfolds
            obs.fold1 <- (n0+1):nobs %% nfolds
            obs.fold[y == 0] <- sample(obs.fold0)
            obs.fold[y == 1] <- sample(obs.fold1)
            obs.fold[obs.fold == 0] <- nfolds
        } else {
            obs.fold <- sample(1:nobs %% nfolds)
            obs.fold[obs.fold == 0] <- nfolds
        }
    } else {
        nfolds <- max(obs.fold)
    }

    cv.args <- list(...)
    # all folds should use the same lambda (kappa) sequence
    cv.args$lambda <- lambda
    cv.args$kappa <- kappa

    if (ncores > 1) {
        doParallel::registerDoParallel(cores = ncores)
        cv.res <- foreach(fold = 1:nfolds, .combine = "rbind", .packages=c("glmtlp")) %dopar% {
            fit.fold <- do.call("glmtlp", c(list(X = X[obs.fold != fold, , drop = FALSE], 
                                                 y = y[obs.fold != fold]), 
                                            cv.args))
            yhat <- predict.glmtlp(fit.fold, X=X[obs.fold == fold, , drop=FALSE], type="response") # ny * nlambda(nkappa)
            loss <- loss.glmtlp(y[obs.fold == fold], yhat, family) # a vector of nlambda (nkappa)
            loss
        }
    } else {
        cv.res <- c()
        for (fold in 1:nfolds) {
            fit.fold <- do.call("glmtlp", c(list(X = X[obs.fold != fold, , drop = FALSE], 
                                                 y = y[obs.fold != fold]), 
                                            cv.args))
            yhat <- predict.glmtlp(fit.fold, X=X[obs.fold == fold, , drop = FALSE], type="response")
            loss <- loss.glmtlp(y[obs.fold == fold], yhat, family)
            cv.res <- rbind(cv.res, loss)
        }
    }
    cv.mean <- apply(cv.res, 2, mean)
    cv.se <- apply(cv.res, 2, sd)
    idx.min <- which.min(cv.mean)[1]  # if there are multiple, output the larger one

    out <- structure(list(call = cv.call, 
                          fit = fit, 
                          obs.fold = obs.fold, 
                          cv.mean = cv.mean, 
                          cv.se = cv.se, 
                          idx.min = idx.min,
                          null.dev = loss.glmtlp(y, rep(mean(y), nobs), family)
                          ), 
                     class = "cv.glmtlp")
    if (penalty == "l0") {
        out$lambda <- lambda
        out$kappa <- kappa
        out$kappa.min <- kappa[idx.min]
    } else {
        out$lambda <- lambda
        out$lambda.min <- lambda[idx.min]
    }
    out
}
