
get.setting.comp <- function(n, p, rho, s, beta.type, snr, seed, log.file, nrep = 5) {
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
          " has been finished.\n",
          sep = ""
        ),
        file = log.file, append = TRUE
        )
        seed <- seed + 1
      }
      params <- matrix(rep(c(n, p, rho, s, beta.type, snr), nrow(comp.res)),
        nrow = nrow(comp.res), byrow = TRUE
      )
      colnames(params) <- c("n", "p", "rho", "s", "beta.type", "snr")
      cbind(params, comp.res)
    },
    error = function(cond) {
      cat(paste("\nThe setting n=", n, ", p=", p, ", s=", s,
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

gen.data <- function(n, p, rho = 0, s = 5, beta.type = 1, snr = 1, seed = 2021) {
  set.seed(seed)

  X <- matrix(rnorm(n * p), n, p)

  if (rho != 0) {
    for (j in 2:p) {
      X[, j] <- sqrt(1 - rho^2) * X[, j] + rho * X[, j - 1]
    }
  }

  # L <- cholesky.ar1.root(rho, p)
  # if (rho != 0) {
  #   X = X %*% t(L)
  # }

  s <- min(s, p)
  beta <- rep(0, p)
  if (beta.type == 1) {
    beta[round(seq(1, p, length = s))] <- 1
  } else if (beta.type == 2) {
    beta[1:s] <- 1
  } else if (beta.type == 3) {
    beta[1:s] <- seq(10, 0.5, length = s)
  } else if (beta.type == 4) {
    beta[1:6] <- c(-10, -6, -2, 2, 6, 10)
  } else {
    beta[1:s] <- 1
    if (s + 1 <= p) {
      beta[(s + 1):p] <- 0.5^(1:(p - s))
    }
  }

  vmu <- sum((  crossprod(cholesky.ar1.root(rho, p), beta) )^2)

  # vmu <- as.numeric(t(beta) %*% L %*% t(L) %*% beta)
  # rm(L)
  sigma <- sqrt(vmu / snr)

  mu <- plogis(as.numeric(X %*% beta + rnorm(n) * sigma))
  y <- rbinom(n, 1, mu)

  list(X = X, y = y)
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

comp.models <- function(data, log.file, cutoff = 50) {
  res.glmnet <- fit.glmnet(data, log.file, cutoff)
  res.ncvreg <- fit.ncvreg(data, log.file, cutoff)
  res.logistic <- fit.logistic(data, log.file, cutoff)

  res <- t(sapply(
    list(res.glmnet, res.ncvreg, res.logistic),
    function(res) {
      if ("measures" %in% names(res)) {
        return(res$measures)
      }
    }
  ))

  methods <- c("glmnet", "ncvreg", "logistic")
  res <- cbind(methods, res)
  res
}

fit.glmnet <- function(data, log.file, cutoff) {
  require(glmnet)
  res <- tryCatch(
    {
      X <- data$X
      y <- data$y
      t_start <- Sys.time()
      mod <- glmnet(
        x = X, y = y, family = "binomial", standardize = FALSE,
        lambda = get.glmnet.lambda(X, y)
      )
      t_end <- Sys.time()
      runtime <- t_end - t_start

      intercept <- as.numeric(mod$a0[cutoff])
      coefs <- as.numeric(mod$beta[, cutoff])
      lambda <- mod$lambda[cutoff]

      measures <- c(runtime, get.obj.val(X, y, intercept, coefs, lambda))
      names(measures) <- c("runtime", "obj")

      list(model = mod, measures = measures)
    },
    error = function(cond) {
      cat(paste("\nThere is something wrong with glmnet.\n", sep = ""),
        file = log.file, append = TRUE
      )
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
      return(NA)
    },
    finally = {
    }
  )
  return(res)
}

get.glmnet.lambda <- function(x, y, nlambda = 100) {
  nobs <- as.integer(nrow(x))
  nvars <- as.integer(ncol(x))
  lambda.min.ratio <- ifelse(nobs <= nvars, 0.01, 1e-04)
  lambda.max <- max(abs(crossprod(x, (y - mean(y))))) / nobs
  lambda <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), len = nlambda))
  lambda
}

get.obj.val <- function(X, y, intercept, coefs, lambda) {
  eta <- intercept + X %*% coefs
  mu.hat <- plogis(eta)
  -mean(y * log(mu.hat) + (1 - y) * log(1 - mu.hat)) + lambda * sum(abs(coefs))
}


fit.ncvreg <- function(data, log.file, cutoff) {
  require(ncvreg)
  res <- tryCatch(
    {
      X <- data$X
      y <- data$y
      t_start <- Sys.time()
      mod <- ncvreg(
        X = X, y = y, family = "binomial",
        lambda = get.glmnet.lambda(X, y), penalty = "lasso",
        returnX = FALSE
      )
      t_end <- Sys.time()
      runtime <- t_end - t_start

      intercept <- as.numeric(mod$beta[1, cutoff])
      coefs <- as.numeric(mod$beta[-1, cutoff])
      lambda <- mod$lambda[cutoff]

      measures <- c(runtime, get.obj.val(X, y, intercept, coefs, lambda))
      names(measures) <- c("runtime", "obj")

      list(model = mod, measures = measures)
    },
    error = function(cond) {
      cat(paste("\nThere is something wrong with ncvreg\n", sep = ""),
        file = log.file, append = TRUE
      )
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
      return(NA)
    },
    finally = {
    }
  )
  return(res)
}

fit.logistic <- function(data, log.file, cutoff) {
  require(glmtlp)
  res <- tryCatch(
    {
      X <- data$X
      y <- data$y
      t_start <- Sys.time()
      mod <- logistic(
        X = X, y = y, lambda = get.glmnet.lambda(X, y),
        delta = 2.0, tol = 1e-4
      )
      t_end <- Sys.time()
      runtime <- t_end - t_start

      intercept <- as.numeric(mod$b0[cutoff])
      coefs <- as.numeric(mod$b[, cutoff])
      lambda <- mod$lambda[cutoff]

      measures <- c(runtime, get.obj.val(X, y, intercept, coefs, lambda))
      names(measures) <- c("runtime", "obj")

      list(model = mod, measures = measures)
    },
    error = function(cond) {
      cat(paste("\nThere is something wrong with logistic in glmtlp.\n", sep = ""),
        file = log.file, append = TRUE
      )
      cat("Here's the original error message:", file = log.file, append = TRUE)
      cat(cond, file = log.file, append = TRUE)
      return(NA)
    },
    finally = {
    }
  )
  return(res)
}