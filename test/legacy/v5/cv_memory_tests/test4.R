# need to restart the session after running one test
# or check Activity Monitor for the memory occupation

library(pryr)
library(glmtlp)  # use the cran published version 2.0.0
library(peakRAM)
library(doParallel)
test_dir <- "/Users/yuyang/Documents/Projects/ongoing/glmtlp-dev/test/v5/cv_memory_tests"
log.file <- file.path(test_dir, "test_log.txt")
cat(paste("\n[", Sys.time(), "] \tPackage version: ", packageVersion("glmtlp"), 
          sep = ""), file = log.file, append = TRUE)

# test 4
mem <- peakRAM(
  {
    # data preparation
    X <- matrix(rnorm(10000 * 1000), 10000, 1000)
    y <- rnorm(10000)
    nfolds <- 10
    family <- "gaussian"
    penalty <- "l1"
    seed <- 2021
    obs.fold <- NULL
    ncores <- 2
    fit <- glmtlp(X = X, y = y, family = family, penalty = penalty) # will do the error checking and generate lambda, kappa sequences
    nobs <- nrow(X)
    family <- fit$family
    penalty <- fit$penalty
    lambda <- fit$lambda
    kappa <- fit$kappa
    if (!is.null(seed)) set.seed(seed)
    if (nfolds > nobs) stop(paste("nfolds (", nfolds, ") cannot be larger than the number of observations (", nobs, ")", sep = ""))
    # generate fold
    if (is.null(obs.fold)) {
      if (family == "binomial") { # stratified sampling for binary data
        n0 <- sum(y == 0)
        obs.fold[y == 0] <- sample(rep(1:nfolds, length.out = n0))
        obs.fold[y == 1] <- sample(rep(1:nfolds, length.out = nobs - n0))
      } else {
        obs.fold <- sample(rep(1:nfolds, length.out = nobs))
      }
    } else {
      nfolds <- max(obs.fold)
    }
    doParallel::registerDoParallel(cores = ncores)
    cv.res <- foreach(fold = 1:nfolds, .combine = "rbind", .packages = c("glmtlp")) %dopar% {
      fit.fold <- glmtlp(X, y, weights = 1 * (obs.fold != fold), 
                         lambda = lambda, kappa = kappa, family = family, penalty = penalty)
    }
  }
)
cat(paste("\nTest 4 Runtime: ", mem$Elapsed_Time_sec, ", \t Peak Memory: ", 
          mem$Peak_RAM_Used_MiB, "\n", sep = ""), file = log.file, append = TRUE)
rm(list = ls())




