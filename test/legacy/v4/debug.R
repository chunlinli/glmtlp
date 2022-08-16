setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("testutils4.R")

#================================= Toy testing ===============================
data <- gen.data(n = 500, p = 1000, rho = 0, s = 5, beta.type = 1, 
                 snr = 6, seed = 2021)

library(glmtlp)
mod <- glmtlp(data$X, data$y, family = "gaussian", penalty = "tlp")
coef(mod)

X <- data$X
y <- data$y
xdim <- dim(X)
nobs <- as.integer(xdim[1])
nvars <- as.integer(xdim[2])
family = "gaussian"
penalty = "tlp"
nlambda=100
lambda.min.ratio=ifelse(nobs < nvars, 1e-2, 1e-4)
lambda=NULL
kappa=NULL
tau=0.3*sqrt(log(nvars)/nobs)
delta=2.0
tol=1e-4
weights=NULL 
penalty.factor=rep(1.0, nvars)
standardize=TRUE
dc.maxit=20
cd.maxit=10000
nr.maxit=500





# check X
xdim <- dim(X)
nobs <- as.integer(xdim[1])
nvars <- as.integer(xdim[2])

if (is.null(xdim) | (xdim[2] <= 1)) stop("X should be a matrix with 2 or more columns")
if (!inherits(X, "matrix")) {
  tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
  if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix")
}
if (typeof(X)=="character") stop("X must be a numeric matrix")
if (typeof(X)=="integer") storage.mode(X) <- "double"



# check y
y <- drop(y) # Delete the dimensions of an array which has only one level.
dimy <- dim(y)
nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
if (nrowy != nobs) stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of X (", nobs, ")", sep = ""))
if (!is.double(y)) {
  op <- options(warn=2)
  on.exit(options(op))
  y <- tryCatch(
    error = function(cond) stop("y must be numeric or able to be coerced to numeric"),
    as.double(y))
  options(op)
}
if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to the model")
if (family == "binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
if (family == "binomial" & !identical(sort(unique(y)), 0:1)) {
  y <- as.double(y == max(y))
} 

# check penalty.factor, weights, and s
penalty.factor <- as.double(penalty.factor) 
if (length(penalty.factor) != nvars) stop(paste("the length of penalty.factor (", length(penalty.factor), ") not equal to the number of variables (", nvars, ")", sep = ""))
if (is.null(weights)) {
  weights <- rep(1.0, nobs)
} else if (length(weights) != nobs) {
  stop(paste("number of elements in weights (", length(weights), ") not equal to the number of rows of X (", nobs, ")", sep = ""))
} else {
  weights <- as.double(weights)
}

## Deprication support
dots <- list(...)
#   if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

## standardize X
if (standardize) {
  X <- std(X)
}

## check/setup lambda and kappa
if (is.null(lambda)) {
  nlambda <- as.integer(nlambda)
  if (nlambda < 2) stop("nlambda must be at least 2")
  if(lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
  lambda <- setup_lambda(X, y, lambda.min.ratio, nlambda)
  user.lambda <- FALSE
} else {
  nlambda <- length(lambda)
  if (nlambda < 2) stop("the length of lambda must be at least 2")
  if (any(lambda < 0)) stop("lambdas should be non-negative")
  if (!is.double(lambda)) lambda <- as.double(lambda)
  lambda <- sort(lambda, decreasing = TRUE)
  user.lambda <- TRUE
}

if (penalty == "l0") {
  if(is.null(kappa)) {
    kappa <- 1:min(nvars, as.integer(nvars/log(nobs)))
    nkappa <- length(kappa)
    user.kappa <- FALSE
  } else {
    if (any(kappa < 0)) stop("kappa should be non-negative")
    if (!is.integer(kappa)) {
      message("Coercing kappa to integers.")
      kappa <- as.integer(kappa) 
    }
    kappa <- sort(kappa, decreasing = FALSE)
    nkappa <- length(kappa)
    user.kappa <- TRUE
  }
}

## check tau, delta, tol, and maxiters: may add more on tau, delta, and tol checking
tau <- as.double(tau)
delta <- as.double(delta)
tol <- as.double(tol)
dc.maxit <- as.integer(dc.maxit)
cd.maxit <- as.integer(cd.maxit)
nr.maxit <- as.integer(nr.maxit)

## fit
if (family == "gaussian") {
  fit <- switch(penalty,
                "l0" = gaussian_l0_reg(X, y, kappa, lambda, tau, weights, penalty.factor, delta, tol, dc.maxit, cd.maxit), 
                "l1" = gaussian_l1_reg(X, y, lambda, weights, penalty.factor, delta, tol, cd.maxit), 
                "tlp" = gaussian_tlp_reg(X, y, lambda, tau, weights, penalty.factor, delta, tol, dc.maxit, cd.maxit)
  ) 
} else if (family == "binomial") {
  fit <- switch(penalty,
                "l0" = logistic_l0_reg(), 
                "l1" = logistic_l1_reg(X, y, lambda, penalty.factor, delta, tol, nr.maxit, cd.maxit), 
                "tlp" = logistic_tlp_reg()
  )
}

## Unstandardize
# beta <- matrix(0, nrow = nvars + 1, ncol = ifelse(penalty == "l0", nkappa, nlambda))
# b0 <- fit$b0
# b <- fit$b
# bb <- b / attr(X, "scaled:scale")
# beta[2:(nvars+1), ] <- bb
# beta[1, ] <- b0 - crossprod(attr(X, "scaled:center"), bb)

intercept <- fit$b0
beta <- fit$b

if (standardize) {
  beta <- beta / attr(X, "scaled:scale")
  intercept <- intercept - crossprod(attr(X, "scaled:center"), beta)
} 


## Names
varnames <- colnames(X)
if (is.null(varnames)) varnames <- paste("V", seq(nvars), sep = "") 
rownames(beta) <- varnames
if (penalty == "l0") {
  colnames(beta) <- paste(kappa)
} else {
  colnames(beta) <- lambda_names(lambda)
}
## Output
out <- structure(list(beta = beta, 
                      call = this.call,
                      family = family,
                      intercept = intercept, 
                      penalty = penalty, 
                      penalty.factor = penalty.factor, 
                      weights = weights
), 
class = "glmtlp")
if (penalty == "l0") {
  out$lambda <- lambda
  out$kappa <- kappa
  out$user.lambda <- user.lambda
  out$user.kappa <- user.kappa
} else if (penalty == "l1") {
  out$lambda <- lambda
  out$user.lambda <- user.lambda
} else {
  out$tau <- tau
}
out