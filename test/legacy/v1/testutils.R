require(glmnet)
require(ncvreg)
require(glmtlp) # the built one

gen.data <- function(n, p, p0, snr, seed=8053) {
  # p0: number of non-zero coefficients
  # snr: signal-to-noise ratio
  set.seed(seed)
  X <- matrix(rnorm(n*p),n,p)
  Z <- rnorm(n)
  for(j in 1:p) {
    X[,j] <- X[,j] + Z
    X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
  }
  set.seed(seed + 1)
  beta <- rbinom(p0, 1, 0.5) * 2 - 1 # 1 or -1
  
  set.seed(seed + 2)
  noise <- rnorm(n)
  f.x <- 1 + 0.5 * X[, 1:p0] %*% beta 
  k <- c(sqrt(var(f.x) / snr / var(noise)))
  noise <- k * noise
  y <- f.x + noise
  
  list(X=X, y=y)
}


comp.models <- function(data, cutoff=50) {
  X <- data$X
  y <- data$y
  fitted <- fit.models(X, y)
  
  obj.glm <- get.obj.glm(X, y, fitted$glm, cutoff)
  obj.ncv <- get.obj.ncv(X, y, fitted$ncv, cutoff)
  obj.lasso <- get.obj.lasso(X, y, fitted$lasso, cutoff)
  
  obj.vals <- c(obj.glm, obj.ncv, obj.lasso, 
                obj.lasso - obj.glm, obj.lasso - obj.ncv)
  names(obj.vals) <- c("glmnet", "ncvreg", "lasso", "lasso-glmnet", "lasso-ncvreg")
  
  list(obj.vals=obj.vals, runtime=fitted$runtime)
}


fit.models <- function(X, y) {
  
  res.glmnet <- fit.glmnet(X, y)
  res.ncvreg <- fit.ncvreg(X, y)
  res.lasso <- fit.lasso(X, y)
  
  t.cost <- c(res.glmnet$runtime, res.ncvreg$runtime, res.lasso$runtime)
  names(t.cost) <- c("glmnet", "ncvreg", "lasso")
  
  list(glm=res.glmnet$model, ncv=res.ncvreg$model, lasso=res.lasso$model, 
       runtime=t.cost)
}

fit.glmnet <- function(X, y) {
  res <- tryCatch(
    {
      t_start <- Sys.time()
      m <- glmnet(x=X, y=y, standardize=FALSE, lambda=get.glmnet.lambda(X, y), thresh=1e-8)
      t_end <- Sys.time()
      t <- t_end - t_start
      list(model=m, runtime=t)
    },
    error=function(cond) {
      message(paste("There is something wrong with glmnet.", sep = ""))
      message("Here's the original error message:")
      message(cond)
      return(list(model=NA, runtime=NA))
    },
    finally={
    }
  )
  return(res)
}

fit.ncvreg <- function(X, y) {
  res <- tryCatch(
    {
      t_start <- Sys.time()
      m <- ncvreg(X=X, y=y, lambda = get.glmnet.lambda(X, y), 
                  penalty = "lasso", returnX = FALSE)
      t_end <- Sys.time()
      t <- t_end - t_start
      list(model=m, runtime=t)
    },
    error=function(cond) {
      message(paste("There is something wrong with ncvreg", sep = ""))
      message("Here's the original error message:")
      message(cond)
      return(list(model=NA, runtime=NA))
    },
    finally={
    }
  )
  return(res)
}

get.glmnet.lambda <- function(x, y, nlambda=100) {
  nobs <- as.integer(nrow(x))
  nvars <- as.integer(ncol(x))
  lambda.min.ratio <- ifelse(nobs <= nvars, 0.01, 1e-04)
  lambda.max <- max(abs(crossprod(x, (y - mean(y))))) / nobs
  lambda <- exp(seq(log(lambda.max),log(lambda.min.ratio*lambda.max),len=nlambda))
  lambda
}

fit.lasso <- function(X, y) {
  res <- tryCatch(
    {
      t_start <- Sys.time()
      m <- lasso(X=X, y=y, lambda=get.glmnet.lambda(X, y), delta=2.0, tol = 1e-4) 
      t_end <- Sys.time()
      t <- t_end - t_start
      list(model=m, runtime=t)
    },
    error=function(cond) {
      message(paste("There is something wrong with ncvreg", sep = ""))
      message("Here's the original error message:")
      message(cond)
      return(list(model=NA, runtime=NA))
    },
    finally={
    }
  )
  return(res)
}

get.obj.glm <- function(X, y, m, cutoff=50) {
  res <- tryCatch(
    {
      a <- as.numeric(m$a0[cutoff])
      b <- as.numeric(m$beta[, cutoff])
      lambda <- m$lambda[cutoff]
      obj.value(X, y, a, b, lambda)
    },
    error=function(cond) {
      message(paste("There is something wrong with glmnet objective calculation.", sep = ""))
      message("Here's the original error message:")
      message(cond)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

obj.value <- function(X, y, intercept, coefs, lambda) {
  mean((y - X%*%coefs - intercept)^2/2) + lambda*sum(abs(coefs))
}

get.obj.ncv <- function(X, y, m, cutoff=50) {
  res <- tryCatch(
    {
      a <- m$beta[1,cutoff]
      b <- as.numeric(m$beta[-1,cutoff])
      lambda <- m$lambda[cutoff]
      obj.value(X, y, a, b, lambda)
    },
    error=function(cond) {
      message(paste("There is something wrong with ncvreg objective calculation.", sep = ""))
      message("Here's the original error message:")
      message(cond)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

get.obj.lasso <- function(X, y, m, cutoff=50) {
  res <- tryCatch(
    {
      a <- m$b0[cutoff]
      b <- m$b[, cutoff]
      lambda <- m$lambda[cutoff]
      obj.value(X, y, a, b, lambda)
    },
    error=function(cond) {
      message(paste("There is something wrong with lasso objective calculation.", sep = ""))
      message("Here's the original error message:")
      message(cond)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}




