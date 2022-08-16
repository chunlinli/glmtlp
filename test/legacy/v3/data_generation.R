
gen_data1 <- function(n, p, rho = 0) {
    X <- matrix(rnorm(n * p), n, p)

    if (rho != 0) {
        for (j in 2:p) {
            X[, j] <- sqrt(1 - rho^2) * X[, j] + rho * X[, j - 1]
        }
    }

    X
}

gen_data2 <- function(n, p, rho = 0) {
    X <- matrix(rnorm(n * p), n, p)

    L <- cholesky.ar1.root(rho, p)
    if (rho != 0) {
        X <- X %*% t(L)
    }

    X
}

cholesky.ar1.root <- function(rho, p) {
    # reference 1: https://blogs.sas.com/content/iml/2018/10/03/ar1-cholesky-root-simulation.html
    # reference 2: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4455603/ (Direct formulation to Cholesky decomposition of a general nonsingular correlation matrix)
    if (rho != 0) {
        L <- matrix(0, nrow = p, ncol = p)
        L[, 1] <- rho^(0:(p - 1))
        c <- sqrt(1 - rho^2)
        cL <- c * L[, 1]
        for (i in 2:p) {
            L[i:p, i] <- cL[1:(p - i + 1)]
        }
    } else {
        L <- diag(1, p)
    }

    L
}

n <- 1000
p <- 10000

t1 <- Sys.time()
X <- gen_data1(n=n,p=p,rho=0.9)
t2 <- Sys.time()
t2-t1

# t1 <- Sys.time()
# Z <- gen_data2(n=n,p=p,rho=0.9)
# t2 <- Sys.time()
# t2-t1


t1 <- Sys.time()
Y <- cholesky.ar1.root(rho=0.9, p=p)
t2 <- Sys.time()
t2 - t1