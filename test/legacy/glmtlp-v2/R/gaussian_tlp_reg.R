#dyn.load("./src/glmtlp")

# remove the initial guess of b.
gaussian_tlp_reg <- function(X, y, lambda, tau, weights, penalty.factor, delta, tol, dc.maxit, cd.maxit) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    
    r <- y - mean(y)
    nlambda <- as.integer(length(lambda))

    b0 <- rep(mean(y), nlambda)
    b <- matrix(0, p, nlambda)

    # call regularized version
    .Call("gaussian_tlp", b0, b, r, X, weights, penalty.factor, lambda, nlambda, tau, 
          n, p, delta, tol, dc.maxit, cd.maxit)

    list(b = b, b0 = b0, lambda = lambda, tau = tau)
}