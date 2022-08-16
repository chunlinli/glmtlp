# remove the initial guess of b.

gaussian_l1_reg <- function(X, y, lambda, weights, penalty.factor, delta, tol, cd.maxit) {
    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))

    r <- y - mean(y)
    nlambda <- as.integer(length(lambda))

    b0 <- rep(mean(y), nlambda)
    b <- matrix(0, p, nlambda)
    
    # call C interface: rlasso.cc
    .Call("gaussian_l1", b0, b, r, X, weights, penalty.factor, lambda, nlambda,
        n, p, delta, tol, cd.maxit)

    list(b = b, b0 = b0, lambda = lambda)
}

# b0 should be a vec: length(b0) == length(lambda)