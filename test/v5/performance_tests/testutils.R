
get.setting.comp <- function(n, p, rho, kappa, beta.type, snr, seed, log.file, nrep = 5) {
  out <- tryCatch(
    {
      comp.res <- c()
      for (i in 1:nrep) {
        data <- gen.gaussian.data(n, p, rho, kappa, beta.type, snr, seed)
        comp.res <- rbind(comp.res, comp.gaussian.models(data, log.file))
        rm(data)
        cat(paste("\n[", Sys.time(), "] The setting n=", n, ", p=", p, 
                  ", kappa=", kappa, ", rho=", rho, ", beta.type=", beta.type, 
                  ", snr=", snr, ", seed=", seed, ", iter=", i, 
                  " has been finished.\n",sep = ""), 
            file = log.file, append = TRUE
        )
        seed <- seed + 1
      }
      params <- matrix(rep(c(n, p, rho, kappa, beta.type, snr), nrow(comp.res)),
        nrow = nrow(comp.res), byrow = TRUE
      )
      colnames(params) <- c("n", "p", "rho", "kappa", "beta.type", "snr")
      cbind(params, comp.res)
    },
    error = function(cond) {
      cat(paste("\nThe setting n=", n, ", p=", p, ", kappa=", kappa,
        ", rho=", rho, ", beta.type=", beta.type,
        ", snr=", snr, ", seed=", seed,
        " has something wrong.\n",
        sep = ""
      ),
      file = log.file, append = TRUE
      )
      cat("Here's the original error message:",
        file = log.file, append = TRUE
      )
      cat(cond,
        file = log.file, append = TRUE
      )
      return(NA)
    },
    finally = {
    }
  )
  return(out)
}


gen.gaussian.data <- function(n, p, rho=0, kappa=5, beta.type=1, snr=1, seed=2021) {
  set.seed(seed)
  X <- matrix(rnorm(n*p), n, p)
  if (rho != 0) {
    for (j in 2:p) {
      X[, j] <- sqrt(1 - rho^2) * X[, j] + rho * X[, j - 1]
    }
  }
  beta <- gen.beta(kappa, p, beta.type)
  vmu <- sum((crossprod(cholesky.ar1.root(rho, p), beta))^2)
  sigma <- sqrt(vmu / snr)
  y <- as.numeric(X %*% beta + rnorm(n) * sigma)
  
  if (p > 1) {
    colnames(X) <- paste("V", seq(p), sep = "") 
  }
  
  list(X = X, y = y, beta = beta, sigma = sigma, rho=rho)
}

gen.beta <- function(kappa, p, beta.type) {
  kappa <- min(kappa, p)
  beta <- rep(0, p)
  if (beta.type == 1) {
    beta[round(seq(1, p, length = kappa))] <- 1
  } else if (beta.type == 2) {
    beta[1:kappa] <- 1
  } else if (beta.type == 3) {
    beta[1:kappa] <- seq(10, 0.5, length = kappa)
  } else if (beta.type == 4) {
    beta[1:6] <- c(-10, -6, -2, 2, 6, 10)
  } else {
    beta[1:kappa] <- 1
    if (kappa + 1 <= p) {
      beta[(kappa + 1):p] <- 0.5^(1:(p - kappa))
    }
  }
  beta
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

comp.gaussian.models <- function(data, log.file) {
  
  res.l0reg <- fit.l0reg(data, log.file)
  res.l1reg <- fit.l1reg(data, log.file)
  res.tlpreg <- fit.tlpreg(data, log.file)
  res.scad <- fit.scad(data, log.file)
  res.mcp <- fit.mcp(data, log.file)
  res.glmnet <- fit.glmnet(data, log.file)
  
  res <- t(sapply(list(res.l0reg, res.l1reg, res.tlpreg, 
                       res.scad, res.mcp, res.glmnet), 
                  function(res) {if ("measures" %in% names(res)) return(res$measures)}
  ))
  
  methods <- c("l0reg", "l1reg", "tlpreg", "scad", "mcp", "glmnet")
  res <- cbind(methods, res)
  res
}

fit.l0reg <- function(data, log.file) {
  require(glmtlp)
  require(parallel)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      cv.fit <- cv.glmtlp(X = data$X, y = data$y, seed = 2021, 
                          family = "gaussian", penalty = "l0", 
                          ncores=min(detectCores() / 2, 10))
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      coefs <- coef(cv.fit)
      beta.hat <- coefs[-1]
      intercept <- coefs[1]
      
      measures <- c(runtime, get.measures(data, beta.hat, intercept))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=cv.fit, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with l0reg.\n", sep = ""), 
          file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(paste(cond), file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

get.measures <- function(data, beta.hat, intercept) {
  # reference: http://arxiv.org/abs/1707.08692 Section 3.1
  p <- ncol(data$X)
  yhat <- data$X %*% beta.hat + intercept
  err <- mean((yhat - data$y)^2)
  delta <- beta.hat - data$beta
  risk <- sum((crossprod(cholesky.ar1.root(data$rho, p), delta))^2) + intercept^2
  err.test <- risk + data$sigma^2
  risk.null <- sum((crossprod(cholesky.ar1.root(rho, p), data$beta))^2)
  err.null <- risk.null + data$sigma^2
  risk.rel <- risk / risk.null
  err.rel <- err.test / data$sigma^2
  prop <- 1 - err.test / err.null
  nzs <- sum(beta.hat != 0)
  tp <- sum((beta.hat != 0) * (data$beta != 0))
  fp <- sum((beta.hat != 0) * (data$beta == 0))
  fn <- sum((beta.hat == 0) * (data$beta != 0))
  F1 <- 2 * tp / (2 * tp + fp + fn)
  hd <- fp + fn
  hd.rel <- hd / sum(data$beta != 0)
  opt <- (err.test - err) / err
  
  measures <- c(err, err.test, prop, risk, nzs, tp, fp, fn, F1, opt, 
                risk.rel, err.rel, hd, hd.rel)
  measures
}


fit.l1reg <- function(data, log.file) {
  require(glmtlp)
  require(parallel)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      cv.fit <- cv.glmtlp(X = data$X, y = data$y, seed = 2021, 
                          family = "gaussian", penalty = "l1", 
                          ncores=min(detectCores() / 2, 10))
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      coefs <- coef(cv.fit)
      beta.hat <- coefs[-1]
      intercept <- coefs[1]
      
      measures <- c(runtime, get.measures(data, beta.hat, intercept))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=cv.fit, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with l1reg.\n", sep = ""), 
          file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(paste(cond), file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

fit.tlpreg <- function(data, log.file) {
  require(glmtlp)
  require(parallel)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      cv.fit <- cv.glmtlp(X = data$X, y = data$y, seed = 2021, 
                          family = "gaussian", penalty = "tlp", 
                          ncores=min(detectCores() / 2, 10))
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      coefs <- coef(cv.fit)
      beta.hat <- coefs[-1]
      intercept <- coefs[1]
      
      measures <- c(runtime, get.measures(data, beta.hat, intercept))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=cv.fit, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with tlpreg.\n", sep = ""), 
          file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(paste(cond), file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

fit.scad <- function(data, log.file) {
  require(ncvreg)
  require(parallel)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      cl <- makeCluster(min(detectCores() / 2, 10))
      mod <- cv.ncvreg(X=data$X, y=data$y, penalty="SCAD", cluster = cl)
      stopCluster(cl)
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      beta.hat <- as.numeric(mod$fit$beta[-1,mod$min])
      intercept <- as.numeric(mod$fit$beta[1,mod$min])
      
      measures <- c(runtime, get.measures(data, beta.hat, intercept))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with scad.\n", sep = ""), 
          file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(paste(cond), file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

fit.mcp <- function(data, log.file) {
  require(ncvreg)
  require(parallel)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      cl <- makeCluster(min(detectCores() / 2, 10))
      mod <- cv.ncvreg(X=data$X, y=data$y, penalty="MCP", cluster = cl)
      stopCluster(cl)
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      beta.hat <- as.numeric(mod$fit$beta[-1,mod$min])
      intercept <- as.numeric(mod$fit$beta[1,mod$min])
      
      measures <- c(runtime, get.measures(data, beta.hat, intercept))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with mcp.\n", sep = ""), 
          file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(paste(cond), file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

fit.glmnet <- function(data, log.file) {
  require(glmnet)
  require(doMC)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      registerDoMC(cores = min(detectCores() / 2, 10))
      mod <- cv.glmnet(x=data$X, y=data$y, alpha=1, parallel = TRUE)
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      beta.hat <- as.numeric(coef(mod, s = "lambda.min"))[-1]
      intercept <- as.numeric(coef(mod, s = "lambda.min"))[1]
      
      measures <- c(runtime, get.measures(data, beta.hat, intercept))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with lasso.\n", sep = ""), 
          file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(paste(cond), file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}


gen.binomial.data <- function(n, p, rho = 0, kappa = 5, beta.type = 1, seed = 2021) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  if (rho != 0) {
    for (j in 2:p) {
      X[, j] <- sqrt(1 - rho^2) * X[, j] + rho * X[, j - 1]
    }
  }
  beta <- gen.beta(kappa, p, beta.type)
  mu <- plogis(as.numeric(X %*% beta))
  y <- rbinom(n, 1, mu)
  
  if (p > 1) {
    colnames(X) <- paste("V", seq(p), sep = "") 
  }
  
  list(X = X, y = y, beta = beta)
  
}