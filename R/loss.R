# /**********
#     R Interface: Loss functions for generalized linear models
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



# need to revise to allow for customed weights
loss_glmtlp <- function(y, yhat, weights = rep(1, length(y)), family) {
    n <- ifelse(is.null(dim(yhat)), length(yhat), dim(yhat)[1])
    if (n != length(y)) {
        stop("dim 1 of yhat must be equal to the length of y")
    }
    yhat <- matrix(yhat, nrow = n)

    if (family == "gaussian") {
        dev <- (y - yhat)^2
    } else if (family == "binomial") {
        dev <- matrix(NA, nrow = nrow(yhat), ncol = ncol(yhat))
        yhat[yhat < 0.00001] <- 0.00001
        yhat[yhat > 0.99999] <- 0.99999
        dev[y == 0, ] <- -2.0 * log(1 - yhat[y == 0, , drop = FALSE])
        dev[y == 1, ] <- -2.0 * log(yhat[y == 1, , drop = FALSE])
    } else if (family == "poisson") {
        dev <- 2.0 * (y * log(y / yhat + 0.00001) - (y - yhat))
    } else {
        stop("family should be one of gaussian, binomial, poisson")
    }
    dev <- apply(dev, 2, sum)
    drop(dev)
}
