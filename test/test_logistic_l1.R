

library(glmnet)
library(ncvreg)
#source('test/test_fns.R')
library(glmtlp)


# DATA GENERATION
#set.seed(2021)


n <- 1000
p <- 10000
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}

mu <- plogis(1 + 0.5*(X[,1] - X[,2] + X[,10] - X[,50] + X[,200]))
y <- rbinom(n,1,mu)

# SPEED TEST: single solution path (100 lambdas)
# Compare with glmnet and ncvreg
t0 <- Sys.time()
m3 <- glmtlp(X=X, y=y, family="binomial", method = "l1")
t1 <- Sys.time()
print("ours")
print(t1 - t0)


t0 <- Sys.time()
m1 <- glmnet(x=X, y=y, lambda = m3$lambda, family = "binomial", standardize=FALSE)
t1 <- Sys.time()
print("glmnet")
print(t1 - t0)


t0 <- Sys.time()
m2 <- ncvreg(X=X, y=y, family="binomial", lambda = m3$lambda, penalty = "lasso",
             returnX=FALSE)
t1 <- Sys.time()
print("ncvreg")
print(t1 - t0)


k = 50
# compute (penalized negative) log-likelihood
lambda <- m1$lambda[k]

b1 <- as.numeric(m1$beta[,k])
a1 <- as.numeric(m1$a0[k])
eta <- a1 + X %*% b1
mu.hat <- plogis(eta)
lik1 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b1))

b2 <- as.numeric(m2$beta[-1,50])
a2 <- as.numeric(m2$beta[1,50])
eta <- a2 + X %*% b2
mu.hat <- plogis(eta)
lik2 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b2))

b3 <- as.numeric(m3$beta[,k])
a3 <- as.numeric(m3$intercept[k])
eta <- a3 + X %*% b3
mu.hat <- plogis(eta)
lik3 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b3))

# the smaller the better
print("likelihood")
print(c(lik1,lik2,lik3))



