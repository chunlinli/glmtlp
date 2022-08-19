# /**********
#     R Interface: Constrained likelihood ratio test for glmtlp
#
#     Copyright (C) 2022 Chunlin Li
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


test_glmtlp <- function(X, y, test_null, penalty_factor = NULL,
                        family = c("gaussian", "binomial", "poisson"),
                        model_selection = c("bic", "cv"), ...) {
    this_call <- match.call(expand.dots = TRUE)
    family <- match.arg(family)
    model_selection <- match.arg(model_selection)

    xdim <- dim(X)
    if (is.null(xdim) || xdim[2] < 1) {
        stop("X should be a matrix")
    }

    nobs <- xdim[1]
    nvars <- xdim[2]
    if (length(test_null) < 1) stop("please provide null hypothesis")
    test_null <- as.integer(unique(test_null))
    if (any(test_null < 1 | test_null > nvars)) {
        stop("tested variable must be between [1, nvars]")
    }
    if (is.null(penalty_factor)) {
        penalty_factor <- rep(1, nvars)
    } else {
        if (any(penalty_factor[test_null] != 0)) {
            message("penalty_factor of tested variable is not 0,
            setting to 0")
        }
    }
    penalty_factor[test_null] <- 0

    if (model_selection == "cv") {
        cv_object <- cv_glmtlp(
            X = X, y = y, family = family,
            penalty_factor = penalty_factor,
            method = "tlp-constrained", ...
        )
        idx_min <- cv_object$idx_min
        dev1 <- cv_object$fit$deviance[idx_min]
        kappa_min <- cv_object$kappa_min
        lambda <- cv_object$lambda
    } else {
        fit1 <- glmtlp(
            X = X, y = y, family = family,
            penalty_factor = penalty_factor,
            method = "tlp-constrained", ...
        )
        bic <- fit1$deviance + log(nrow(X)) * fit1$kappa
        idx_min <- which.min(bic)[1]
        dev1 <- fit1$deviance[idx_min]
        kappa_min <- fit1$kappa[idx_min]
        lambda <- fit1$lambda
    }
    fit0 <- glmtlp(
        X = X, y = y, family = family,
        penalty_factor = penalty_factor,
        lambda = lambda,
        kappa = kappa_min,
        method = "tlp-constrained",
        exclude = test_null,
        ...
    )
    dev0 <- fit0$deviance[idx_min]

    # test statistic for gaussian is not defined in this way
    log_like_ratio <- ifelse(family == "gaussian",
        (dev0 - dev1) * (nobs - kappa_min - 1) / dev1,
        dev0 - dev1
    )
    pvalue <- pchisq(log_like_ratio,
        df = length(test_null),
        lower.tail = FALSE
    )

    structure(list(
        call = this_call,
        pvalue = pvalue,
        log_like_ratio = log_like_ratio,
        dev1 = dev1,
        dev0 = dev0,
        df = length(test_null),
        kappa_min = kappa_min,
        test_null = test_null,
        model_selection = model_selection,
        penalty_factor = penalty_factor
    ),
    class = "test_glmtlp"
    )
}
