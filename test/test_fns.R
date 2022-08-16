library(Rcpp)
library(RcppEigen)
sourceCpp("src/glm.cpp")
sourceCpp("src/sum.cpp")

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
        fit <- glm_solver(
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


setup_lambda <- function(X, y, weights, lambda_min_ratio, nlambda) {
    lambda_max <- get_lambda_max(X, y, weights)
    lambda <- exp(seq(
        from = log(lambda_max),
        to = log(lambda_min_ratio * lambda_max),
        length.out = nlambda
    ))
    lambda
}

get_lambda_max <- function(X, y, weights) {
    rw <- (y - weighted.mean(y, weights)) * weights
    max(abs(crossprod(X, rw)), na.rm = TRUE) / nrow(X)
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

lambda_names <- function(l) {
    if (length(l) > 1) {
        d <- ceiling(-log10(-max(diff(l))))
        d <- min(max(d, 4), 10)
    } else {
        d <- 4
    }
    formatC(l, format = "f", digits = d)
}



sumtlp <- function(XX, Xy, nobs,
                   method = c("tlp-constrained", "tlp-regularized", "l1-regularized"),
                   nlambda = ifelse(method == "tlp-constrained", 50, 100),
                   lambda_min_ratio = ifelse(nobs < nvars, 5e-2, 1e-3),
                   lambda = NULL, kappa = NULL,
                   df_max = min(nvars, 40),
                   tau = 0.5 * sqrt(log(nvars) / nobs),
                   delta = 2.0, tol = 1e-4,
                   penalty_factor = rep(1.0, nvars),
                   cd_maxit = 10000, ...) {

    # Coersion
    this_call <- match.call()
    method <- match.arg(method)

    nobs <- as.integer(nobs)
    nvars <- as.integer(length(Xy))

    # check data XX, Xy
    xxdim <- dim(XX)
    if (is.null(xxdim) || any(xxdim != nvars)) {
        stop("X should be a squared matrix")
    }

    if (any(is.na(Xy)) || any(is.na(XX))) {
        stop("Missing data (NA's) detected. Take actions
                (e.g., removing cases, removing features, imputation)
                to eliminate missing data before passing
                XX and Xy to the model")
    }

    # check penalty_factor
    if (method == "tlp-constrained") penalty_factor <- (penalty_factor != 0) * 1
    penalty_factor <- as.double(penalty_factor)
    if (length(penalty_factor) != nvars) {
        stop(paste("the length of penalty_factor (",
            length(penalty_factor), ") not equal to the number of variables (",
            nvars, ")",
            sep = ""
        ))
    }

    ## check/setup lambda and kappa
    if (is.null(lambda)) {
        nlambda <- as.integer(nlambda)
        if (nlambda < 2) stop("nlambda must be at least 2")
        if (lambda_min_ratio >= 1) {
            stop("lambda_min_ratio should be less than 1")
        }
        lambda <- setup_lambda_sum(Xy, lambda_min_ratio, nlambda)
    } else {
        nlambda <- as.integer(length(lambda))
        if (nlambda < 1) stop("the length of input lambda must be at least 1")
        if (!is.double(lambda)) lambda <- as.double(lambda)
        lambda <- sort(lambda, decreasing = TRUE)
        if (lambda[nlambda] < 0) stop("lambdas should be non-negative")
        lambda_max <- max(abs(Xy))
        if (lambda[1] < lambda_max) {
            lambda <- c(lambda_max, lambda)
        } else {
            if (lambda[nlambda] >= lambda_max) {
                message("null model, try smaller lambdas")
                return(get_null_output_sum(
                    this_call,
                    family, method, penalty_factor
                ))
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

    fit <- sum_solver(
        XX, Xy, penalty_factor, kappa, lambda, delta, tau, tol, cd_maxit, method
    )

    varnames <- colnames(XX)
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
            method = method,
            penalty_factor = penalty_factor,
            loss = fit$loss,
            call = this_call
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
    out
}



cv_glmtlp <- function(X, y, ...,
                      nfolds = 10, obs_fold = NULL,
                      seed = NULL, ncores = 1) {
    cv_call <- match.call(expand.dots = TRUE)

    fit <- glmtlp(X = X, y = y, ...)
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
                family = family, method = method
            )
            yhat <- predict.glmtlp(fit_fold,
                X = X[obs_fold == fold, , drop = FALSE],
                type = "response"
            )
            loss <- loss_glmtlp(y[obs_fold == fold], yhat, family)
            loss
        }
    } else {
        cv_res <- c()
        for (fold in 1:nfolds) {
            fit_fold <- glmtlp(X, y,
                weights = 1 * (obs_fold != fold),
                lambda = lambda, kappa = kappa,
                family = family, method = method
            )
            yhat <- predict.glmtlp(fit_fold,
                X = X[obs_fold == fold, , drop = FALSE],
                type = "response"
            )
            loss <- loss_glmtlp(y[obs_fold == fold], yhat, family)
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
        null_dev = loss_glmtlp(y, rep(mean(y), nobs), family)
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



# H0 is index to test: null is "beta[H0] == 0"
# this function needs improvement:
# it is not efficient
# it is misses corner cases

test_glmtlp <- function(X, y, H0, penalty_factor = NULL,
                        model_selection = c("bic", "cv"), ...) {
    this_call <- match.call(expand.dots = TRUE)
    model_selection <- match.arg(model_selection)

    xdim <- dim(X)
    if (is.null(xdim) || xdim[2] < 1) {
        stop("X should be a matrix")
    }

    nvars <- xdim[2]
    if (length(H0) < 1) stop("please provide null hypothesis H0")
    H0 <- as.integer(unique(H0))
    if (any(H0 < 1 | max(H0) > nvars)) stop("H0 must be between 1:nvars")
    if (is.null(penalty_factor)) {
        penalty_factor <- rep(1, nvars)
    } else {
        if (any(penalty_factor[H0] != 0)) {
            message("penalty_factor of tested var is not 0, setting it to 0")
        }
    }
    penalty_factor[H0] <- 0

    if (model_selection == "cv") {
        cv_object <- cv_glmtlp(
            X = X, y = y,
            penalty_factor = penalty_factor,
            method = "tlp-constrained", ...
        )
        idx_min <- cv_object$idx_min
        dev1 <- cv_object$fit$deviance[idx_min]
        kappa_min <- cv_object$kappa_min
        lambda <- cv_object$lambda
    } else {
        fit1 <- glmtlp(
            X = X, y = y,
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
        X = X[, -H0], y = y,
        penalty_factor = penalty_factor[-H0],
        lambda = lambda,
        kappa = kappa_min,
        method = "tlp-constrained", ...
    )
    dev0 <- fit0$deviance[idx_min]
    test_stat <- dev0 - dev1
    pvalue <- pchisq(test_stat, df = length(H0), lower.tail = FALSE)

    structure(list(
        call = this_call,
        pvalue = pvalue,
        test_stat = test_stat,
        dev1 = dev1,
        dev0 = dev0,
        kappa_min = kappa_min,
        H0 = H0,
        model_selection = model_selection,
        penalty_factor = penalty_factor
    ),
    class = "test_glmtlp"
    )
}
