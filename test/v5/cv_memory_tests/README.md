# Memory Consumption Test

This test targets at checking the memory consumption by different cv code writing methods.

`test_log.txt` includes the memory test results.
- `2.0.1` is the new version updated on Nov 24 and `2.0.0` is the CRAN version.
- on both 2.0.0 and 2.0.1, test 4 occupies the least RAM.
- on 2.0.1, test 4 occupies about half the memory as on 2.0.0, while the other tests don't change much.
- all 4 tests have been run 3 times on both versions to reduce the effect of randomness.
- the difference among the 4 tests lies in how cv is performed. Check the details below.



## Test 1
```r
# test 1
cv.args <- list(penalty=penalty, family=family)
cv.args$lambda <- lambda
cv.args$kappa <- kappa
cv.res <- foreach(fold = 1:nfolds, .combine = "rbind", .packages = c("glmtlp")) %dopar% {
    fit.fold <- do.call("glmtlp", c(
    list(X = X[obs.fold != fold, , drop = FALSE], y = y[obs.fold != fold]), cv.args
    ))
}
```

## Test 2
```r
# test 2
cv.args <- list(penalty=penalty, family=family)
cv.args$lambda <- lambda
cv.args$kappa <- kappa
cv.res <- foreach(fold = 1:nfolds, .combine = "rbind", .packages = c("glmtlp")) %dopar% {
    cv.args$X <- X[obs.fold != fold, , drop = FALSE]
    cv.args$y <- y[obs.fold != fold]
    fit.fold <- do.call("glmtlp", cv.args)
}
```

## Test 3
```r
# test 3
cv.args <- list(penalty=penalty, family=family)
cv.args$lambda <- lambda
cv.args$kappa <- kappa
cv.args$X <- X
cv.args$y <- y
cv.res <- foreach(fold = 1:nfolds, .combine = "rbind", .packages = c("glmtlp")) %dopar% {
    cv.args$weights <- 1 * (obs.fold != fold)
    fit.fold <- do.call("glmtlp", cv.args)
}
```

## Test 4
```r
# test 4
cv.res <- foreach(fold = 1:nfolds, .combine = "rbind", .packages = c("glmtlp")) %dopar% {
      fit.fold <- glmtlp(X, y, weights = 1 * (obs.fold != fold), 
                         lambda = lambda, kappa = kappa, family = family, penalty = penalty)
    }
```