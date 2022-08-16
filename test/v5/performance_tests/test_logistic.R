

library(glmnet)
library(ncvreg)
library(glmtlp)

# DATA GENERATION
set.seed(2021)
n <- 500
p <- 400
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}

mu <- plogis(1 + 0.5*(X[,1] - X[,2] + X[,10] - X[,50] + X[,200]))  
y <- rbinom(n,1,mu)


runtime <- c()
# SPEED TEST: single solution path (100 lambdas) 
# Compare with glmnet and ncvreg
t0 <- Sys.time()
m1 <- glmnet(x=X, y=y, family = "binomial", standardize=FALSE)
t1 <- Sys.time()
runtime <- c(runtime, difftime(t1, t0, units = 'secs'))


t0 <- Sys.time()
m2 <- ncvreg(X=X, y=y, family="binomial", lambda = m1$lambda, penalty = "lasso", 
             returnX=FALSE)
t1 <- Sys.time()
runtime <- c(runtime, difftime(t1, t0, units = 'secs'))

t0 <- Sys.time()
m3 <- glmtlp(X=X, y=y, tol = 1e-4, family="binomial", penalty="l1") 
t1 <- Sys.time()
runtime <- c(runtime, difftime(t1, t0, units = 'secs'))

t0 <- Sys.time()
m4 <- glmtlp(X=X, y=y, tol = 1e-4, family="binomial", penalty="l0") 
t1 <- Sys.time()
runtime <- c(runtime, difftime(t1, t0, units = 'secs'))

t0 <- Sys.time()
m5 <- glmtlp(X=X, y=y, tol = 1e-4, family="binomial", penalty="tlp") 
t1 <- Sys.time()
runtime <- c(runtime, difftime(t1, t0, units = 'secs'))

# compute (penalized negative) log-likelihood
lambda <- m1$lambda[50]

b1 <- as.numeric(m1$beta[,50])
a1 <- as.numeric(m1$a0[50])
eta <- a1 + X %*% b1
mu.hat <- plogis(eta)
lik1 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b1))

b2 <- as.numeric(m2$beta[-1,50])
a2 <- as.numeric(m2$beta[1,50])
eta <- a2 + X %*% b2
mu.hat <- plogis(eta)
lik2 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b2))

b3 <- as.numeric(m3$beta[,50])
a3 <- as.numeric(m3$intercept[50])
eta <- a3 + X %*% b3
mu.hat <- plogis(eta)
lik3 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b3))

b4 <- as.numeric(m4$beta[,50])
a4 <- as.numeric(m4$intercept[50])
eta <- a4 + X %*% b4
mu.hat <- plogis(eta)
lik4 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b4))

b5 <- as.numeric(m5$beta[,50])
a5 <- as.numeric(m5$intercept[50])
eta <- a5 + X %*% b5
mu.hat <- plogis(eta)
lik5 <- -mean(y * log(mu.hat) + (1-y)*log(1-mu.hat)) + lambda*sum(abs(b5))

# the smaller the better
lik <- c(lik1, lik2, lik3, lik4, lik5)
names(lik) <- c("glmnet", "ncvreg", "l1", "l0", "tlp")
names(runtime) <- c("glmnet", "ncvreg", "l1", "l0", "tlp")
lik
runtime

m5
