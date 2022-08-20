

library(glmtlp)

set.seed(1110)
n <- 1000
p <- 5000
X <- matrix(rnorm(n * p), n, p)
Z <- rnorm(n)
for (j in 1:p) {
    X[, j] <- X[, j] + Z
    X[, j] <- (X[, j] - mean(X[, j])) / sd(X[, j])
}
y <- 1 + 0.5 * (X[, 1] - X[, 2] + X[, 10] - X[, 50] + X[, 200]) + rnorm(n)


m <- test_glmtlp(X, y,
    test_null = c(3, 11),
    family = "gaussian", method = "tlp-c",
    ncores = 4, model_selection = "cv"
)

m$pvalue