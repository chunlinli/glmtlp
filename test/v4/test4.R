setwd(".")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./testutils4.R")

#================================= Toy testing ===============================
data <- gen.data(n = 500, p = 1000, rho = 0, s = 5, beta.type = 1, 
                 snr = 6, seed = 2021)
library(glmtlp)
fit <- glmtlp(data$X, data$y, family = "binomial", penalty = "l1")
coef(fit)
predict(fit, X=data$X)
plot(fit, xvar = "log_lambda", label = FALSE)

fit <- glmtlp(data$X, data$y, family = "gaussian", penalty = "l0")
coef(fit)
predict(fit, X=data$X)
plot(fit, xvar = "kappa")
plot(fit, xvar = "log_kappa")

fit <- glmtlp(data$X, data$y, family = "gaussian", penalty = "l1")
coef(fit)
predict(fit, X=data$X)
plot(fit, xvar = "log_lambda", label = FALSE)

fit <- glmtlp(data$X, data$y, family = "gaussian", penalty = "tlp")
coef(fit)
predict(fit, X=data$X)
plot(fit, xvar = "log_lambda", label = FALSE)



t_start <- Sys.time()
cv.fit <- cv.glmtlp(data$X, data$y, family = "binomial", penalty = "l1", ncores=5)
t_end <- Sys.time()
runtime <- t_end - t_start
cv.fit
runtime
plot(cv.fit)

t_start <- Sys.time()
cv.fit <- cv.glmtlp(data$X, data$y, family = "gaussian", penalty = "l0", ncores=5)
t_end <- Sys.time()
runtime <- t_end - t_start
cv.fit
runtime
plot(cv.fit)

t_start <- Sys.time()
cv.fit <- cv.glmtlp(data$X, data$y, family = "gaussian", penalty = "l1", ncores=5)
t_end <- Sys.time()
runtime <- t_end - t_start
cv.fit
runtime

t_start <- Sys.time()
cv.fit <- cv.glmtlp(data$X, data$y, family = "gaussian", penalty = "tlp", ncores=5)
t_end <- Sys.time()
runtime <- t_end - t_start
cv.fit
runtime



