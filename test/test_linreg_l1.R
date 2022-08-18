

library(glmnet)
library(ncvreg)
#source('test/test_fns.R')
library(glmtlp)

# DATA GENERATION
set.seed(1110)
n <- 500
p <- 10000
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- 1 + 0.2*(X[,1] - X[,2] + X[,10] - X[,50] + X[,200]) + rnorm(n)


# SPEED TEST: single solution path (100 lambdas)
# Compare with glmnet and ncvreg
t0 <- Sys.time()
m3 <- glmtlp(X=X, y=y, method="l1", family = "gaussian", ncores = 4)
t1 <- Sys.time()
print("ours")
print(t1-t0)

t0 <- Sys.time()
m1 <- glmnet(x=X, y=y, family="gaussian", lambda = m3$lambda, standardize=FALSE)
t1 <- Sys.time()
print("glmnet")
print(t1-t0)

t0 <- Sys.time()
m2 <- ncvreg(X=X, y=y, lambda = m3$lambda, penalty = "lasso", returnX = FALSE)
t1 <- Sys.time()
print("ncvreg")
print(t1-t0)

# ACCURACY TEST: for example, let's compare the 50th solution
a1 <- as.numeric(m1$a0[50])
a2 <- m2$beta[1,50]
a3 <- m3$intercept[50]
b1 <- as.numeric(m1$beta[, 50])
b2 <- as.numeric(m2$beta[-1,50])
b3 <- m3$beta[, 50]

lambda <- m1$lambda[50]

# OBJECTIVE VALUE: the smaller the better
print(mean((y - X%*%b3 - a3)^2/2) + lambda*sum(abs(b3)))
print(mean((y - X%*%b1 - a1)^2/2) + lambda*sum(abs(b1)))
print(mean((y - X%*%b2 - a2)^2/2) + lambda*sum(abs(b2)))


