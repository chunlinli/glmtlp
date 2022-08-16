
source('test/test_fns.R')

# DATA GENERATION
set.seed(1110)
n <- 500
p <- 100
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- 1 + 0.3*(X[,1] - X[,2] + X[,10] - X[,50]) + rnorm(n)
y <- y - mean(y)

# for sumtlp, XX must be correlation matrix
# Xy must satisfy that X is standardized (mean 0, sd 1) and y is centered (mean 0).
XX <- cov(X)
Xy <- crossprod(X,y)/n

# TEST: single solution path (100 lambdas)
t0 <- Sys.time()
m1 <- glmtlp(X=X, y=y, method="tlp-r", family = "g")
t1 <- Sys.time()
print("glmtlp: computing 100 lambdas")
print(t1-t0)

t0 <- Sys.time()
m2 <- sumtlp(XX=XX, Xy=Xy, nobs = n, method="tlp-r")
t1 <- Sys.time()
print("sumtlp: computing 100 lambdas")
print(t1-t0)

k <- 31
idx1 <- which(m1$beta[, k] != 0)
idx2 <- which(m2$beta[, k] != 0)

b1 <- m1$beta[idx1,k]
b2 <- m2$beta[idx2,k]

print("glmtlp selects:")
print(b1)
print("sumtlp selects:")
print(b2)

t0 <- Sys.time()
m3 <- glmtlp(X=X, y=y, method="tlp-c", family = "g", kappa = 1:20)
t1 <- Sys.time()
print("glmtlp: computing 100 lambdas")
print(t1-t0)

t0 <- Sys.time()
m4 <- sumtlp(XX=XX, Xy=Xy, nobs = n, method="tlp-c",kappa=1:20)
t1 <- Sys.time()
print("sumtlp: computing 100 lambdas")
print(t1-t0)

k <- 5
idx <- which(abs(m3$beta[,k])!=0)
b <- m3$beta[idx,k]
print("glmtlp selects:")
print(b)

idx <- which(abs(m4$beta[,k])!=0)
b <- m4$beta[idx,k]
print("sumtlp selects:")
print(b)




