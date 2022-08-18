#' Predict Method for a "cv.glmtlp" Object.
#'
#' @description
#' Makes predictions for a cross-validated glmtlp model, using
#'   the stored \code{"glmtlp"} object, and the optimal value chosen for
#'   \code{lambda}.
#'
#' @aliases coef.cv.glmtlp predict.cv.glmtlp
#' @param object Fitted \code{"cv.glmtlp"} object.
#' @param X X Matrix of new values for \code{X} at which predictions are to be
#'   made. Must be a matrix.
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
#'   sequence at which predictions are to be made. Default is the \code{idx.min}
#'   stored in the \code{cv.glmtp} object.
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
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' cv.fit <- cv.glmtlp(X, y, family = "gaussian", penalty = "l1")
#' predict(cv.fit, X = X[1:5, ])
#' coef(cv.fit)
#' predict(cv.fit, X = X[1:5, ], lambda = 0.1)
#'
#' @method predict cv.glmtlp
#' @export
#' @export predict.cv.glmtlp

predict.cv_glmtlp <- function(object, X, type=c("link","response","class","coefficients","numnzs","varnzs"),
                              lambda=NULL, kappa=NULL, which=object$idx_min, ...) {
  type <- match.arg(type)
  predict.glmtlp(object$fit, X=X, type=type, lambda=lambda, kappa=kappa,
                 which=which, ...)
}

#' Extract coefficients from a cv_glmtlp object
#'
#' @method coef cv_glmtlp
#' @rdname predict.cv_glmtlp
#' @export
#' @export coef.cv_glmtlp
#'
coef.cv_glmtlp <- function(object, lambda=NULL, kappa=NULL, which=object$idx_min, ...) {
  coef.glmtlp(object$fit, lambda=lambda, kappa=kappa, which=which, ...)
}
