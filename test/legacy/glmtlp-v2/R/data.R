#' A simulated binomial dataset.
#'
#' A dataset simulated for illustrating logistic regression models.
#'
#' @format A list with two elements: design matrix \code{X} and \code{y}.
#' \describe{
#'   \item{X}{design matrix}
#'   \item{y}{response}
#'   ...
#' }
#' 
#' @usage data(bin_data)
#' 
#' @examples
#' data("bin_data")
#' cv.fit <- cv.glmtlp(bin_data$X, bin_data$y, family = "binomial", penalty = "l1")
#' plot(cv.fit)
#' 
"bin_data"

#' A simulated gaussian dataset.
#'
#' A dataset simulated for illustrating linear regression models.
#'
#' @format A list with five elements: design matrix \code{X}, response \code{y}, 
#'   correlation structure of the covariates \code{Sigma}, true beta \code{beta}, 
#'   and the noise level \code{sigma}.
#' \describe{
#'   \item{X}{design matrix}
#'   \item{y}{response}
#'   \item{Sigma}{correlation matrix of the covariates}
#'   \item{beta}{true beta values}
#'   \item{sigma}{the noise level}
#'   ...
#' }
#' 
#' @usage data(gau_data)
#' 
#' @examples 
#' data("gau_data")
#' cv.fit <- cv.glmtlp(gau_data$X, gau_data$y, family = "gaussian", penalty = "tlp")
#' plot(cv.fit)
#' 
"gau_data"