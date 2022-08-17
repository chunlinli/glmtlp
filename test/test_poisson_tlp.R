
library(tictoc)
library(glmnet)
library(ncvreg)
#source('test/test_fns.R')
library(glmtlp)


# DATA GENERATION
#set.seed(1117)

n <- 2000
p <- 20000
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}

mu <- exp(1 + .4*(X[,1] - X[,2] + X[,10] - X[,50] + X[,200]))
y <- rpois(n,mu)

# SPEED TEST: single solution path (100 lambdas)
# Compare with glmnet and ncvreg
tic("TLP-Regularized")
m1 <- glmtlp(X=X, y=y, family = "poisson", method="tlp-r")
toc()

k <- 50
beta <- m1$beta[,k]
idx <- which(beta!=0)
b1 <- beta[idx]

print("TLP-regularized selects:")
print(b1)


tic("TLP-Constrained")
m2 <- glmtlp(X=X, y=y, family="poisson", method="tlp-c")
toc()

k <- 6
idx <- which(abs(m2$beta[,k])!=0)
b2 <- m2$beta[idx,k]
print("TLP-constrained selects:")
print(b2)





# tic("MCP")
# m3 <- ncvreg(X=X, y=y, family="poisson", penalty="MCP",returnX=FALSE)
# toc()
#
# k <- 32
# beta <- m3$beta[-1,k]
# idx <- which(beta != 0)
# b3 <- beta[idx]
#
# print("MCP selects:")
# print(b3)
#
#
# print(b2/b3)

