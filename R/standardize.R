#' @importFrom stats model.matrix
std <- function(X) {
  if (typeof(X) == 'integer') storage.mode(X) <- 'double'
  if (!inherits(X, "matrix")) {
    if (is.numeric(X)) {
      X <- matrix(as.double(X), ncol=1)
    } else {
      tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
  }
  
  scale <- apply(X, 2, sd_)

  sX <- scale(X, center = TRUE, scale = scale)
  dimnames(sX) <- dimnames(X)

  sX
}

sd_ <- function(y) sqrt(sum((y-mean(y))^2)/length(y)) 
