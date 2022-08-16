
setwd(".")
source("testutils2.R")

get.setting.extra.metrics <- function(n, p, rho, s, beta.type, snr, seed, log.file, nrep=5) {
  out <- tryCatch(
    {
      extra.res <- c()
      for (i in 1:nrep) {
        data <- gen.data(n, p, rho, s, beta.type, snr, seed)
        extra.res <- rbind(extra.res, get.extra.two.metrics(data))
        rm(data)
        cat(paste("\n[", Sys.time(), "] The setting n=", n, ", p=", p, ", s=", s,
                  ", rho=", rho, ", beta.type=", beta.type, 
                  ", snr=", snr, ", seed=", seed, ", iter=", i, 
                  " has been finished.\n", sep = ""),
            file = log.file, append = TRUE)
        seed <- seed + 1
      }
      params <- matrix(rep(c(n, p, rho, s, beta.type, snr), nrow(extra.res)), 
                       nrow = nrow(extra.res), byrow = TRUE)
      colnames(params) <- c("n", 'p', 'rho', 's', 'beta.type', 'snr')
      cbind(params, extra.res)
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

get.extra.two.metrics <- function(data) {
  
  sigma2 <- data$sigma^2
  risk.null <- c(data$beta %*% data$Sigma %*% data$beta)
  res <- matrix(rep(c(sigma2, risk.null), 5), nrow=5, byrow = T)
  
  methods <- c("l0reg", "tlpreg", "scad", "mcp", "lasso")
  res <- cbind(methods, res)
  res
}

require(doSNOW)

n.candidates <- c(100, 200, 500, 1000)
p.candidates <- c(10, 100, 1000)
rho.candidates <- c(0, 0.35, 0.7)
s.candidates <- c(5, 10, 20)
beta.type.candidates <- c(1:3, 5)
snr.candidates <- exp(seq(log(0.05),log(6),length=10))
seed <- 2021
nrep <- 5
result.file <- "results/extra_metrics.csv"
log.file <- "results/log_extra_metrics.txt"
if (file.exists(log.file)) { file.remove(log.file) }
if (file.exists(result.file)) { file.remove(result.file) }

cat(paste(paste(rep("=", 80), collapse = ""), "\n\n", 
          "Parameters:", "\n\n", 
          "n.candidates = ", paste(n.candidates, collapse = ", "), ".\n\n", 
          "p.candidates = ", paste(p.candidates, collapse = ", "), ".\n\n", 
          "rho.candidates = ", paste(rho.candidates, collapse = ", "), ".\n\n", 
          "s.candidates = ", paste(s.candidates, collapse = ", "), ".\n\n", 
          "beta.type.candidates = ", paste(beta.type.candidates, collapse = ", "), ".\n\n", 
          "snr.candidates = ", paste(snr.candidates, collapse = ", "), ".\n\n", 
          "seed = ", seed, ".\n\n", 
          "nrep = ", nrep, ".\n\n", 
          "result.file = ", result.file, ".\n\n", 
          "log.file = ", log.file, ".\n\n", 
          paste(rep("=", 80), collapse = ""), "\n",
          sep = ""), 
    file = log.file, append = TRUE)

res <- c() 
cl <- makeCluster(5)
registerDoSNOW(cl)
for (n in n.candidates) {
  for (p in p.candidates) {
    for (s in s.candidates) {
      for (rho in rho.candidates) {
        for (beta.type in beta.type.candidates) {
          comp.res <- foreach(snr = snr.candidates, .combine = rbind, 
                              .packages = c("glmnet", "ncvreg", "glmtlp")) %dopar% {
                                get.setting.extra.metrics(n, p, rho, s, beta.type, snr, seed, log.file, nrep)
                              }
          seed <- seed + nrep + 1
          res <- rbind(res, comp.res)
          colnames(res) <- c("n", 'p', 'rho', 's', 'beta.type', 'snr',
                             "method", "sigma2", "risk.null")
          write.csv(res, result.file, row.names = FALSE)
          cat(paste("\n[", Sys.time(), "] Update test results to ", result.file, ".\n\n", sep = ""), 
              file = log.file, append = TRUE)
        }
      }
    }
  }
}
stopCluster(cl)





df1 <- read.csv("results/tests02.csv", header = T)
df2 <- read.csv("results/extra_metrics.csv", header = T)
c(nrow(df1), ncol(df1))
c(nrow(df2), ncol(df2))
sum(df1[,1:7] != df2[, 1:7])

df1$risk.rel <- df1$risk / df2$risk.null
df1$err.rel <- df1$err.test / df2$sigma2
write.csv(df1, "results/tests02_appended.csv")