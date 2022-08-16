
# assumes X has been standardized, y has been centered if gaussian

setup_lambda <- function(X, y, lambda.min.ratio, nlambda) {
    n <- as.integer(nrow(X))

    r <- y - mean(y)
    lambda.max <- max(abs(crossprod(X, r)), na.rm = TRUE) / n

    lambda <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), length.out = nlambda))
    lambda
}

