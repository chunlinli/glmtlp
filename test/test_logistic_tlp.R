
library(tictoc)
library(ncvreg)
#library(glmtlp)
source('test/test_fns.R')


# DATA GENERATION
set.seed(1110)
n <- 5000
p <- 4000
X <- matrix(rnorm(n * p), n, p)
Z <- rnorm(n)
for (j in 1:p) {
  X[, j] <- X[, j] + Z
  X[, j] <- (X[, j] - mean(X[, j])) / sd(X[, j])
}

mu <- plogis(1 + 0.5 * (X[, 1] - X[, 2] + X[, 10] - X[, 50] + X[, 200]))
y <- rbinom(n, 1, mu)


tic("TLP-Regularized")
m1 <- glmtlp(X=X, y=y, family = "binomial", method="tlp-r", cd_maxit = 20000)
toc()

k <- 30
beta <- m1$beta[,k]
idx <- which(beta!=0)
b <- beta[idx]

print("TLP-regularized selects:")
print(b)


tic("TLP-Constrained")
m2 <- glmtlp(X=X, y=y, family="binomial", method="tlp-c", cd_maxit = 20000)
toc()

k <- 6
idx <- which(abs(m2$beta[,k])!=0)
b <- m2$beta[idx,k]
print("TLP-constrained selects:")
print(b)


tic("MCP")
m3 <- ncvreg(X=X, y=y, family="binomial", penalty="MCP",returnX=FALSE)
toc()

k <- 32
beta <- m3$beta[-1,k]
idx <- which(beta != 0)
b <- beta[idx]

print("MCP selects:")
print(b)


library(glmnet)
tic("TLP-Regularized")
m <- glmnet(x=X, y=y, family = "binomial", lambda.min.ratio = 0.001)
toc()
