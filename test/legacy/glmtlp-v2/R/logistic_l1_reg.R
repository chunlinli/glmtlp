#dyn.load("./src/glmtlp")

# remove the initial guess of b.

logistic_l1_reg <- function(X, y, lambda, penalty.factor, delta, tol, nr.maxit, cd.maxit) {
    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))

    nlambda <- as.integer(length(lambda))

    b0 <- rep(0, nlambda)
    b <- matrix(0, p, nlambda)
    
    # call C interface: rlasso.cc
    .Call("logistic_l1", b0, b, y, X, penalty.factor, lambda, nlambda,
        n, p, delta, tol, nr.maxit, cd.maxit)

    list(b = b, b0 = b0, lambda = lambda)
}

# b0 should be a vec: length(b0) == length(lambda)
