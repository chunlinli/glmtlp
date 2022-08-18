# /**********
#     R Interface: Generate lambda sequence.
#
#     Copyright (C) 2021-2022 Yu Yang, Chunlin Li
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

#' Generate lambda sequence.
#'
#' @param X Input matrix, of dimension \code{nobs} x \code{nvars};
#'   each row is  an observation vector.
#' @param y Response variable, of length \code{nobs}. For \code{family="gaussian"},
#'   it should be quantitative; for \code{family="binomial"}, it should be either
#'   a factor with two levels or a binary vector.
#' @param weights Observation weights.
#' @param lambda_min_ratio The smallest value for \code{lambda}, as a fraction of
#'   \code{lambda.max}, the smallest value for which all coefficients are zero.
#'   The default depends on the sample size \code{nobs} relative to the number
#'   of variables \code{nvars}.
#' @param nlambda The number of \code{lambda} values.
#'
#' @importFrom stats weighted.mean
#'
setup_lambda <- function(X, y, weights, lambda_min_ratio, nlambda) {
    lambda_max <- get_lambda_max(X, y, weights)
    lambda <- exp(seq(
        from = log(lambda_max),
        to = log(lambda_min_ratio * lambda_max),
        length.out = nlambda
    ))
    lambda
}

# this must be written in C++ for big matrix
get_lambda_max <- function(X, y, weights) {
    rw <- (y - weighted.mean(y, weights)) * weights
    max(abs(crossprod(X, rw))) / nrow(X)
}

setup_lambda_sum <- function(Xy, lambda_min_ratio, nlambda) {
  lambda_max <- max(abs(Xy))
  lambda <- exp(seq(
    from = log(lambda_max),
    to = log(lambda_min_ratio * lambda_max),
    length.out = nlambda
  ))
  lambda
}