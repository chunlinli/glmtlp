#' Predict Method for a "glmtlp" Object
#'
#' @description
#' Predicts fitted values, logits, coefficients and more from a fitted
#'   \code{glmtlp} object.
#'
#' @details
#' \code{coef(...)} is equivalent to \code{predict(type="coefficients",...)}
#'
#' @aliases coef.glmtlp predict.glmtlp
#'
#' @param object Fitted \code{glmtlp} model object.
#' @param X Matrix of new values for \code{X} at which predictions are to be
#'   made. Must be a matrix. This argument will not used for
#'   \code{type=c("coefficients","numnz", "varnz")}.
#' @param type Type of prediction to be made. For \code{"gaussian"} models, type
#'   \code{"link"} and \code{"response"} are equivalent and both give the fitted
#'   values. For \code{"binomial"} models, type \code{"link"} gives the linear
#'   predictors and type \code{"response"} gives the fitted probabilities.
#'   Type \code{"coefficients"} computes the coefficients at the provided values
#'   of \code{lambda} or \code{kappa}. Note that for \code{"binomial"}
#'   models, results are returned only for the class corresponding to the second
#'   level of the factor response. Type \code{"class"} applies only to
#'   \code{"binomial"} models, and gives the class label corresponding to the
#'   maximum probability. Type \code{"numnz"} gives the total number of non-zero
#'   coefficients for each value of \code{lambda} or \code{kappa}. Type
#'   \code{"varnz"} gives a list of indices of the nonzero coefficients for
#'   each value of \code{lambda} or \code{kappa}.
#' @param lambda Value of the penalty parameter \code{lambda} at which predictions
#'   are to be made Default is NULL.
#' @param kappa Value of the penalty parameter \code{kappa} at which predictions
#'   are to be made. Default is NULL.
#' @param which Index of the penalty parameter \code{lambda} or \code{kappa}
#'   sequence at which predictions are to be made. Default are the indices for the
#'   entire penalty parameter sequence.
#' @param drop Whether or not keep the dimension that is of length 1.
#' @param \dots Additional arguments.
#'
#' @return The object returned depends on \code{type}.
#'
#' @author Chunlin Li, Yu Yang, Chong Wu
#'   \cr Maintainer: Yu Yang \email{yang6367@umn.edu}
#'
#' @seealso \code{print}, \code{predict}, \code{coef} and \code{plot} methods,
#' and the \code{cv.glmtlp} function.
#'
#' @references Shen, X., Pan, W., & Zhu, Y. (2012).
#'   \emph{Likelihood-based selection and sharp parameter estimation.
#'   Journal of the American Statistical Association, 107(497), 223-232.}
#'   \cr Shen, X., Pan, W., Zhu, Y., & Zhou, H. (2013).
#'   \emph{On constrained and regularized high-dimensional regression.
#'   Annals of the Institute of Statistical Mathematics, 65(5), 807-832.}
#'   \cr Li, C., Shen, X., & Pan, W. (2021).
#'   \emph{Inference for a Large Directed Graphical Model with Interventions.
#'   arXiv preprint arXiv:2110.03805.}
#'   \cr Yang, Y., & Zou, H. (2014).
#'   \emph{A coordinate majorization descent algorithm for l1 penalized learning.
#'   Journal of Statistical Computation and Simulation, 84(1), 84-95.}
#'   \cr Two R package Github: \emph{ncvreg} and \emph{glmnet}.
#'
#' @keywords models regression
#'
#' @examples
#'
#' # Gaussian
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' fit <- glmtlp(X, y, family = "gaussian", penalty = "l1")
#' predict(fit, X = X[1:5, ])
#' coef(fit)
#' predict(fit, X = X[1:5, ], lambda = 0.1)
#'
#' # Binomial
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- sample(c(0,1), 100, replace = TRUE)
#' fit <- glmtlp(X, y, family = "binomial", penalty = "l1")
#' coef(fit)
#' predict(fit, X = X[1:5, ], type = "response")
#' predict(fit, X = X[1:5, ], type = "response", lambda = 0.01)
#' predict(fit, X = X[1:5, ], type = "class", lambda = 0.01)
#' predict(fit, X = X[1:5, ], type = "numnz", lambda = 0.01)
#'
#' @importFrom stats plogis
#' @method predict glmtlp
#' @export
#' @export predict.glmtlp

predict.glmtlp <- function(object, X, type=c("link", "response", "class", "coefficients", "numnz", "varnz"),
                           lambda=NULL, kappa=NULL,
                           which=1:(ifelse(object$method == "constrained-tlp", length(object$kappa), length(object$lambda))), ...) {
  type <- match.arg(type)
  coefs <- coef.glmtlp(object, lambda=lambda, kappa=kappa, which=which, drop=FALSE)
  if (type=="coefficients") return(coefs)

  intercept <- coefs[1,]
  beta <- coefs[-1, , drop=FALSE]

  if (type=="numnz") return(apply(beta!=0, 2, sum))
  if (type=="varnz") return(apply(beta!=0, 2, FUN=which))
  eta <- sweep(X %*% beta, 2, intercept, "+")
  if (type=="link" || object$family=="gaussian") return(drop(eta))
  resp <- switch(object$family,
                 binomial = plogis(eta),
                 poisson = exp(eta))
  if (type=="response") return(drop(resp))
  if (type=="class") {
    if (object$family=="binomial") {
      return(drop(1*(eta>0)))
    } else {
      stop("type='class' can only be used with family='binomial'", call.=FALSE)
    }
  }
}


#' Extract coefficients from a glmtlp object
#'
#' @method coef glmtlp
#' @rdname predict.glmtlp
#' @importFrom stats approx
#' @export
#' @export coef.glmtlp
#'
coef.glmtlp <- function(object, lambda=NULL, kappa=NULL,
                        which=1:(ifelse(object$method == "constrained-tlp",
                                        length(object$kappa), length(object$lambda))),
                        drop=TRUE, ...) {
  if (object$method == "constrained-tlp") {
    if (!is.null(kappa)) {
      if (max(kappa) > max(object$kappa) | min(kappa) < min(object$kappa)) {
        stop('Supplied kappa value(s) are outside the range of the model fit.', call.=FALSE)
      }
      ind <- approx(object$kappa, seq(object$kappa), kappa)$y
      l <- floor(ind)
      r <- ceiling(ind)
      w <- ind %% 1
      beta <- (1-w)*object$beta[, l, drop=FALSE] + w*object$beta[, r, drop=FALSE]
      intercept <- (1-w)*object$intercept[l] + w*object$intercept[r]
      colnames(beta) <- kappa
    }
    else {
      beta <- object$beta[, which, drop=FALSE]
      intercept <- object$intercept[which]
    }
  } else {
    if (!is.null(lambda)) {
      if (max(lambda) > max(object$lambda) | min(lambda) < min(object$lambda)) {
        stop('Supplied lambda value(s) are outside the range of the model fit.', call.=FALSE)
      }
      ind <- approx(object$lambda, seq(object$lambda), lambda)$y
      l <- floor(ind)
      r <- ceiling(ind)
      w <- ind %% 1
      beta <- (1-w)*object$beta[, l, drop=FALSE] + w*object$beta[, r, drop=FALSE]
      intercept <- (1-w)*object$intercept[l] + w*object$intercept[r]
      colnames(beta) <- lambda_names(lambda)
    }
    else {
      beta <- object$beta[, which, drop=FALSE]
      intercept <- object$intercept[which]
    }
  }

  if (drop) {
    return (drop(rbind(intercept, beta)))
  } else {
    return (rbind(intercept, beta))
  }
}
