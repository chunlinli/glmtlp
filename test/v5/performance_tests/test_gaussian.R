setwd(".")
source("testutils.R")

require(doSNOW)

n.candidates <- c(100, 200, 500, 1000)
p.candidates <- c(10, 100, 1000)
rho.candidates <- c(0, 0.35, 0.7)
kappa.candidates <- c(5, 10, 20, 50)
beta.type.candidates <- 1:3
snr.candidates <- exp(seq(log(0.05),log(6),length=10))
seed <- 1127
nrep <- 5
result.file <- "results/tests01.csv"
log.file <- "results/log01.txt"
if (file.exists(log.file)) { file.remove(log.file) }
if (file.exists(result.file)) { file.remove(result.file) }

cat(paste(paste(rep("=", 80), collapse = ""), "\n\n", 
          "Parameters:", "\n\n", 
          "n.candidates = ", paste(n.candidates, collapse = ", "), ".\n\n", 
          "p.candidates = ", paste(p.candidates, collapse = ", "), ".\n\n", 
          "rho.candidates = ", paste(rho.candidates, collapse = ", "), ".\n\n", 
          "kappa.candidates = ", paste(kappa.candidates, collapse = ", "), ".\n\n", 
          "beta.type.candidates = ", paste(beta.type.candidates, collapse = ", "), ".\n\n", 
          "snr.candidates = ", paste(snr.candidates, collapse = ", "), ".\n\n", 
          "seed = ", seed, ".\n\n", 
          "nrep = ", nrep, ".\n\n", 
          "package version = ", packageVersion("glmtlp"), ".\n\n", 
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
    for (kappa in kappa.candidates) {
      for (rho in rho.candidates) {
        for (beta.type in beta.type.candidates) {
          comp.res <- foreach(snr = snr.candidates, .combine = rbind, 
             .packages = c("glmnet", "ncvreg", "glmtlp")) %dopar% {
               get.setting.comp(n, p, rho, kappa, beta.type, snr, seed, log.file, nrep)
             }
          seed <- seed + nrep + 1
          res <- rbind(res, comp.res)
          colnames(res) <- c("n", 'p', 'rho', 'kappa', 'beta.type', 'snr', "method", 
                             "runtime", "err", "err.test", "prop", "risk", "nzs", 
                             "tp", "fp", "fn", "F1", "opt", "risk.rel", "err.rel", 
                             "hd", "hd.rel")
          write.csv(res, result.file, row.names = FALSE)
          cat(paste("\n[", Sys.time(), "] Update test results to ", result.file, ".\n\n", sep = ""), 
              file = log.file, append = TRUE)
        }
      }
    }
  }
}
stopCluster(cl)
rm(list = ls())
