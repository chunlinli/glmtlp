
library(tictoc)
library(glmnet)
library(ncvreg)
source('test/test_fns.R')
#library(glmtlp)


# DATA GENERATION
#set.seed(2021)

n <- 2000
p <- 500
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}

mu <- exp(1 + 0.5*(X[,1] - X[,2] + X[,10] - X[,50] + X[,200]))
y <- rpois(n,mu)

# SPEED TEST: single solution path (100 lambdas)
# Compare with glmnet and ncvreg
tic("TLP-Regularized")
m1 <- glmtlp(X=X, y=y, family = "poisson", method="tlp-r")
toc()

k <- 30
beta <- m1$beta[,k]
idx <- which(beta!=0)
b <- beta[idx]

print("TLP-regularized selects:")
print(b)


tic("TLP-Constrained")
m2 <- glmtlp(X=X, y=y, family="poisson", method="tlp-c")
toc()

k <- 6
idx <- which(abs(m2$beta[,k])!=0)
b <- m2$beta[idx,k]
print("TLP-constrained selects:")
print(b)


tic("MCP")
m3 <- ncvreg(X=X, y=y, family="poisson", penalty="MCP",returnX=FALSE)
toc()

k <- 32
beta <- m3$beta[-1,k]
idx <- which(beta != 0)
b <- beta[idx]

print("MCP selects:")
print(b)

