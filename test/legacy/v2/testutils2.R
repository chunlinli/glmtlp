
get.setting.comp <- function(n, p, rho, s, beta.type, snr, seed, log.file, nrep=5) {
  out <- tryCatch(
    {
      comp.res <- c()
      for (i in 1:nrep) {
        data <- gen.data(n, p, rho, s, beta.type, snr, seed)
        comp.res <- rbind(comp.res, comp.models(data, log.file))
        rm(data)
        cat(paste("\n[", Sys.time(), "] The setting n=", n, ", p=", p, ", s=", s,
                ", rho=", rho, ", beta.type=", beta.type, 
                ", snr=", snr, ", seed=", seed, ", iter=", i, 
                " has been finished.\n", sep = ""),
          file = log.file, append = TRUE)
	      seed <- seed + 1
      }
      params <- matrix(rep(c(n, p, rho, s, beta.type, snr), nrow(comp.res)), 
                       nrow = nrow(comp.res), byrow = TRUE)
      colnames(params) <- c("n", 'p', 'rho', 's', 'beta.type', 'snr')
      cbind(params, comp.res)
    },
    error=function(cond) {
      cat(paste("\nThe setting n=", n, ", p=", p, ", s=", s,
                ", rho=", rho, ", beta.type=", beta.type, 
                ", snr=", snr, ", seed=", seed, 
                " has something wrong.\n", sep = ""), 
          file = log.file, append = TRUE)
      cat("Here's the original error message:", 
          file = log.file, append = TRUE)
      cat(cond,
          file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(out)
}

gen.data <- function(n, p, rho=0, s=5, beta.type=1, snr=1, seed=2021) {
  set.seed(seed)
  X = matrix(rnorm(n*p),n,p)

  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    X = X %*% Sigma.half
  } else {
    Sigma = diag(1,p)
  }
  
  s = min(s,p)
  beta = rep(0,p)
  if (beta.type==1) {
    beta[round(seq(1,p,length=s))] = 1
  } else if (beta.type==2) {
    beta[1:s] = 1
  } else if (beta.type==3) {
    beta[1:s] = seq(10,0.5,length=s)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    if (s+1 <= p) {
      beta[(s+1):p] = 0.5^(1:(p-s))
    }
  }
  
  vmu = as.numeric(t(beta) %*% Sigma %*% beta)
  sigma = sqrt(vmu/snr)
  
  y = as.numeric(X %*% beta + rnorm(n)*sigma)
  
  list(X=X,y=y,Sigma=Sigma,beta=beta,sigma=sigma)
}


comp.models <- function(data, log.file) {
  
  res.l0reg <- fit.l0reg(data, log.file)
  res.tlpreg <- fit.tlpreg(data, log.file)
  res.scad <- fit.scad(data, log.file)
  res.mcp <- fit.mcp(data, log.file)
  res.lasso <- fit.lasso(data, log.file)
  
  res <- t(sapply(list(res.l0reg, res.tlpreg, res.scad, res.mcp, res.lasso), 
                  function(res) {if ("measures" %in% names(res)) return(res$measures)}
    ))
  
  methods <- c("l0reg", "tlpreg", "scad", "mcp", "lasso")
  res <- cbind(methods, res)
  res
}

fit.l0reg <- function(data, log.file) {
  require(glmtlp)
  require(parallel)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      mod <- cv.l0reg(X=data$X, y=data$y, ncores=min(detectCores() / 2, 10))
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      beta.hat <- mod$m$b[, mod$s==mod$s.min]
      beta.hat0 <- mod$m$b0[mod$s==mod$s.min]

      measures <- c(runtime, get.measures(data, beta.hat, beta.hat0))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with l0reg.\n", sep = ""), 
        file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

get.measures <- function(data, beta.hat, beta.hat0) {
  
  yhat <- data$X %*% beta.hat + beta.hat0
  err <- mean((yhat - data$y)^2)
  delta <- beta.hat - data$beta
  risk <- c(t(delta) %*% data$Sigma %*% delta) + beta.hat0^2
  err.test <- risk + data$sigma^2
  risk.null <- c(data$beta %*% data$Sigma %*% data$beta)
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


fit.tlpreg <- function(data, log.file) {
  require(glmtlp)
  require(parallel)
  res <- tryCatch(
    {
      t_start <- Sys.time()
      mod <- cv.tlpreg(X=data$X, y=data$y, ncores=min(detectCores() / 2, 10))
      t_end <- Sys.time()
      runtime <- t_end - t_start
      
      beta.hat <- mod$m$b[,mod$lambda==mod$lambda.min]
      beta.hat0 <- mod$m$b0[mod$lambda==mod$lambda.min]
      
      measures <- c(runtime, get.measures(data, beta.hat, beta.hat0))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with tlpreg.\n", sep = ""), 
        file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
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
      beta.hat0 <- as.numeric(mod$fit$beta[1,mod$min])
      
      measures <- c(runtime, get.measures(data, beta.hat, beta.hat0))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with scad.\n", sep = ""), 
        file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
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
      beta.hat0 <- as.numeric(mod$fit$beta[1,mod$min])
      
      measures <- c(runtime, get.measures(data, beta.hat, beta.hat0))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with mcp.\n", sep = ""), 
        file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

fit.lasso <- function(data, log.file) {
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
      beta.hat0 <- as.numeric(coef(mod, s = "lambda.min"))[1]
      
      measures <- c(runtime, get.measures(data, beta.hat, beta.hat0))
      names(measures) <- c("runtime", "err", "err.test", "prop", "risk", "nzs", 
                           "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                           "hd", "hd.rel")
      
      list(model=mod, measures=measures)
    },
    error=function(cond) {
      cat(paste("\nThere is something wrong with lasso.\n", sep = ""), 
        file = log.file, append = TRUE)
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
      return(NA)
    },
    finally={
    }
  )
  return(res)
}

