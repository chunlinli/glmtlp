#' make predictions from a "cv.glmtlp" object.
#'
#' This function makes predictions from a cross-validated glmtlp model, using
#' the stored \code{"glmtlp"} object, and the optimal value chosen for
#' \code{lambda}.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @aliases coef.cv.glmtlp predict.cv.glmtlp
#' @param object Fitted \code{"cv.glmtlp"} object.
#' @param X Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix. See documentation for \code{predict.glmtlp}.
#' @param type Type of prediction required. Type \code{"link"} gives the linear
#'   predictors for \code{"binomial"}
#'   models; for \code{"gaussian"} models it gives the fitted
#'   values. Type \code{"response"} gives the fitted probabilities for
#'   \code{"binomial"}; for \code{"gaussian"} type \code{"response"} is equivalent 
#'   to type \code{"link"}. Type \code{"coefficients"} computes the coefficients 
#'   at the requested values for \code{lambda} or \code{kappa}.  Note that for \code{"binomial"} 
#'   models, results are returned only for the class corresponding to the second 
#'   level of the factor response. Type \code{"class"} applies only to 
#'   \code{"binomial"} models, and produces the class 
#'   label corresponding to the maximum probability. Type \code{"numnz"} returns the 
#'   total number of non-zeros coefficients for each value of \code{lambda} or 
#'   \code{kappa}. Type \code{"varnz"} returns 
#'   a list of the indices of the nonzero coefficients for each value of 
#'   \code{lambda} or \code{kappa}.
#' @param lambda Value of the penalty parameter \code{lambda} at which predictions
#'   are required. Default is NULL.
#' @param kappa Value of the penalty parameter \code{kappa} at which predictions 
#'   are required. Default is NULL.
#' @param which Index of the penalty parameter \code{lambda} or \code{kappa} 
#'   sequence at which predictions are required. Default is the \code{idx.min} stored
#'   on the CV \code{object}. 
#' @param \dots Additional arguments.
#' 
#' @return The object returned depends on the \dots{} argument which is passed
#' on to the \code{predict} method for \code{glmtlp} objects.
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
#'   \cr Yang, Y., & Zou, H. (2014). \emph{A coordinate majorization descent algorithm 
#'   for l1 penalized learning. Journal of Statistical Computation and 
#'   Simulation, 84(1), 84-95.}
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

predict.cv.glmtlp <- function(object, X, type=c("link","response","class","coefficients","numnzs","varnzs"), 
                              lambda=NULL, kappa=NULL, which=object$idx.min, ...) {
  type <- match.arg(type)
  predict.glmtlp(object$fit, X=X, type=type, lambda=lambda, kappa=kappa, 
                 which=which, ...)
}

#' Extract coefficients from a cv.glmtlp object
#'
#' @method coef cv.glmtlp
#' @rdname predict.cv.glmtlp
#' @export
#' @export coef.cv.glmtlp
#' 
coef.cv.glmtlp <- function(object, lambda=NULL, kappa=NULL, which=object$idx.min, ...) {
  coef.glmtlp(object$fit, lambda=lambda, kappa=kappa, which=which, ...)
}