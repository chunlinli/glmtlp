
# deviance is a generic function, so use loss to represent it in the function names
loss.glmtlp <- function(y, yhat, family) {
    n <- ifelse(is.null(dim(yhat)), length(yhat), dim(yhat)[1])
    if (n != length(y)) stop("Dim 1 of yhat should be the same as the length of y.")
    # convert yhat to a matrix (if given a vector) so we can use apply later
    yhat <- matrix(yhat, nrow = n)

    if (family == "gaussian") {
        Dev <- (y - yhat)^2
    } else if (family == "binomial") {
        Dev <- matrix(NA, nrow = nrow(yhat), ncol = ncol(yhat))
        yhat[yhat < 0.00001] <- 0.00001
        yhat[yhat > 0.99999] <- 0.99999
        Dev[y == 0, ] <- -2 * log(1 - yhat[y == 0, , drop = FALSE]) # drop=FALSE keeps the dimension
        Dev[y == 1, ] <- -2 * log(yhat[y == 1, , drop = FALSE]) 
    }
    dev <- apply(Dev, 2, sum)
    drop(dev) # reduce the dimension if length is 1
    dev
}