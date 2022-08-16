setwd(".")
source("testutils.R")

require(doSNOW)

n.candidates <- c(200, 500, 1000, 5000)[1]
p.candidates <- c(100, 1000, 5000, 50000, 200000)[1]
p0.candidates <- c(5, 10, 20, 50)[1]
snr.candidates <- c(0.1, 0.5, 1, 2, 4)[1]
seed <- 8053
result.file <- "results/tests6.csv"
log.file <- "results/log6.txt"
if (file.exists(result.file)) { file.remove(result.file) }
if (file.exists(log.file)) { file.remove(log.file) }

res <- c()
stop <- FALSE
cl <- makeCluster(5)
registerDoSNOW(cl)
for (n in n.candidates) {
  for (p in p.candidates) {
    if(n == 10000 & p == 100000) {
      stop <- TRUE
      break
    }
    for (p0 in p0.candidates) {
      snr.res <- foreach(snr = snr.candidates, .combine = rbind, 
                         .packages = c("glmnet", "ncvreg", "glmtlp")) %dopar% {
         temp.res <- tryCatch(
           {
             data <- gen.data(n, p, p0, snr, seed)
             comp.res <- comp.models(data)
             rm(data)
             seed <- seed + 3
             c(n, p, p0, snr, unlist(comp.res))
           },
           error=function(cond) {
             cat(paste("The setting n=", n, ", p=", p, ", p0=", p0, 
                       ", snr=", snr, 
                       " has something wrong.", sep = ""), 
                 file = log.file, append = TRUE)
             cat("Here's the original error message:", 
                 file = log.file, append = TRUE)
             cat(cond, 
                 file = log.file, append = TRUE)
           },
           finally={
             cat(paste("The setting n=", n, ", p=", p, ", p0=", p0, 
                       ", snr=", snr, 
                       " has been finished.\n", sep = ""),
                 file = log.file, append = TRUE)
           }
         )    
         temp.res
       }
      res <- rbind(res, snr.res)
      colnames(res) <- c("n", "p", "p0", "snr", 
                         "obj.glmnet", "obj.ncvreg", "obj.lasso", 
                         "obj.lasso-obj.glmnet", "obj.lasso-obj.ncvreg", 
                         "runtime.glmnet", "runtime.ncvreg", "runtime.lasso")
      write.csv(res, result.file, row.names = FALSE)
      cat(paste("Update test results to ", result.file, ".\n\n", sep = ""), 
          file = log.file, append = TRUE)
    }
  }
  if (stop) break
}
stopCluster(cl)