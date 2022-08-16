gaussian_l0_reg <- function(X, y, s, lambda, tau, weights, penalty.factor, delta, tol, dc.maxit, cd.maxit) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))

    # initialize K sequence
    ns <- as.integer(length(s))
    
    r <- y - mean(y)
    nlambda <- as.integer(length(lambda))

    b0 <- rep(mean(y), ns)
    b <- matrix(0, p, ns)

    # call regularized version
    .Call("gaussian_l0", b0, b, r, X, weights, penalty.factor, s, ns, lambda, nlambda, tau, 
          n, p, delta, tol, dc.maxit, cd.maxit)

    list(b = b, b0 = b0, s = s, lambda = lambda, tau = tau)
}