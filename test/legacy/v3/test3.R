setwd(".")
source("testutils3.R")

is_parallel <- TRUE
seed <- 2021
nrep <- 3

n.candidates <- c(100, 200, 500, 1000)
p.candidates <- c(100, 1000, 5000, 40000)  # 40000 * 40000 * 8 B = 12.8 GB
rho.candidates <- c(0, 0.35, 0.7, 0.95)
s.candidates <- c(10, 20, 50)
beta.type.candidates <- 1:3
snr.candidates <- exp(seq(log(0.05), log(6), length = 5))

result.file <- "results/tests07.csv"
log.file <- "results/log07.txt"
if (file.exists(log.file)) {
  file.remove(log.file)
}
if (file.exists(result.file)) {
  file.remove(result.file)
}

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
  "R version = ", R.Version()$version.string, ".\n\n",
  paste(rep("=", 80), collapse = ""), "\n",
  sep = ""
),
file = log.file, append = TRUE
)

res <- c()

if (is_parallel) {
  require(doSNOW)
  cl <- makeCluster(10)
  registerDoSNOW(cl)
  for (n in n.candidates) {
    for (p in p.candidates) {
      for (s in s.candidates) {
        for (rho in rho.candidates) {
          for (beta.type in beta.type.candidates) {
            comp.res <- foreach(
              snr = snr.candidates, .combine = rbind,
              .packages = c("glmnet", "ncvreg", "glmtlp")
            ) %dopar% {
              get.setting.comp(n, p, rho, s, beta.type, snr, seed, log.file, nrep)
            }
            seed <- seed + nrep + 1
            res <- rbind(res, comp.res)
            colnames(res) <- c(
              "n", "p", "rho", "s", "beta.type", "snr", "method",
              "runtime", "obj"
            )
            write.csv(res, result.file, row.names = FALSE)
            cat(paste("\n[", Sys.time(), "] Update test results to ", result.file, ".\n\n", sep = ""),
              file = log.file, append = TRUE
            )
          }
        }
      }
    }
  }
  stopCluster(cl)
} else {
  for (n in n.candidates) {
    for (p in p.candidates) {
      for (s in s.candidates) {
        for (rho in rho.candidates) {
          for (beta.type in beta.type.candidates) {
            comp.res <- c()

            for (snr in snr.candidates) {
              comp.res <- rbind(
                comp.res,
                get.setting.comp(n, p, rho, s, beta.type, snr, seed, log.file, nrep)
              )
            }
            seed <- seed + nrep + 1
            res <- rbind(res, comp.res)
            colnames(res) <- c(
              "n", "p", "rho", "s", "beta.type", "snr", "method",
              "runtime", "obj"
            )
            write.csv(res, result.file, row.names = FALSE)
            cat(paste("\n[", Sys.time(), "] Update test results to ", result.file, ".\n\n", sep = ""),
              file = log.file, append = TRUE
            )
          }
        }
      }
    }
  }
}