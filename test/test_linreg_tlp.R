library(tictoc)
library(ncvreg)
library(glmtlp)
#source('test/test_fns.R')

set.seed(1110)
n <- 1000
p <- 500000
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j] - mean(X[,j])) / sd(X[,j])
}
y <- 1 + .4 * (X[,1] - X[,2] + X[,10] - X[,50] + X[,200]) + rnorm(n)


tic("TLP-Regularized")
m1 <- glmtlp(X=X, y=y, family = "g", method="tlp-r", ncores = 4)
toc()

k <- 40
idx <- which(m1$beta[,k]!=0)
b <- m1$beta[idx,k]

print("TLP-regularized selects:")
print(b)

tic("TLP-Constrained")
m2 <- glmtlp(X=X, y=y, family = "g", method="tlp-c", ncores = 4)
toc()

k <- 6
idx <- which(abs(m2$beta[,k])!=0)
b <- m2$beta[idx,k]

print("TLP-constrained selects:")
print(b)

tic("MCP")
m3 <- ncvreg(X=X, y=y, penalty="MCP",returnX=FALSE)
toc()

k <- 40
beta <- m3$beta[-1,k]
idx <- which(beta != 0)
b <- beta[idx]

print("MCP selects:")
print(b)


















