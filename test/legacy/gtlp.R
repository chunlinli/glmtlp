library(glasso)
## This function implements precision matrix estimation via TLP penalty
gtlp <- function(s,lambda,tau,nobs=NULL,zero=NULL,thr=1.0e-4,dc.maxit=20,cd.maxit=1e4,
                 approx=FALSE,penalize.diagonal=FALSE,w.init=NULL,wi.init=NULL){
  p <- nrow(s)
  if(is.null(wi.init)){
    wi.next <- glasso(s=s,rho=lambda/tau,nobs=nobs,zero=zero,thr=thr,maxit=cd.maxit,
                      approx=approx,penalize.diagonal=penalize.diagonal)$wi
    wi.curr <- wi.next
  }else{
    wi.next <- wi.init
    wi.curr <- wi.next
  }
  for (it in 1:dc.maxit){
    rho.mat <- matrix(1,p,p)*(abs(wi.next)<tau)*lambda/tau
    diag(rho.mat) <- rep(0,p)
    wi.next <- glasso(s=s,rho=rho.mat,nobs=nobs,zero=zero,thr=thr,maxit=cd.maxit,
                      approx=approx,penalize.diagonal=penalize.diagonal,wi.init=wi.next)$wi
    if (mean(abs(wi.next-wi.curr))<thr*(sum(abs(s))-sum(diag(abs(s))))/(p*(p-1))) break
    wi.curr <- wi.next
  }
  if (it == dc.maxit) warning("The difference of convex algorithm for gtlp does not converge.")
  return(wi.next)
}

## This function conducts cross-validation for parameters tunning for precision matrix estimation via TLP penalty
cv.gtlp <- function(Z,lambda.list,tau.list,n.fold,nobs=NULL,zero=NULL,thr=1.0e-4,dc.maxit=20,cd.maxit=1e4,
                    approx=FALSE,penalize.diagonal=FALSE,w.init=NULL,wi.init=NULL){
  tau.len <- length(tau.list)
  lambda.len <- length(lambda.list)
  cvm <- matrix(0,tau.len,lambda.len)
  p <- ncol(Z)
  n <- nrow(Z)
  for (i in 1:tau.len){
    tau <- tau.list[i]
    for (j in 1:lambda.len){
      lambda <- lambda.list[j]
      cvm.sub <- numeric(n.fold)
      folds <- cvFolds(n,n.fold)
      for (k in 1:n.fold){
        Z.train <- Z[folds$subsets[folds$which!=k],]
        Z.test <- Z[folds$subsets[folds$which==k],]
        n1 <- nrow(Z.train)
        n2 <- nrow(Z.test)
        s.train <- cov(Z.train)*(n1-1)/n1
        s.test <- cov(Z.test)*(n2-1)/n2
        wi <- gtlp(s.train,lambda,tau)
        cvm.sub[k] <- 0.5*n2*(log(det(wi))-sum(diag(s.test%*%wi)))
      }
      cvm[i,j] <- sum(cvm.sub)
    }
  }
  try(locations <- as.numeric(which(cvm==min(cvm,na.rm = TRUE),arr.ind = TRUE)[1,]),{warning("You may change the penalty parameters used in gtlp");
    locations <- c(ceiling(length(tau.list)/2),ceiling(length(gamma.list)/2))})
  tau <- tau.list[locations[1]]
  lambda <- lambda.list[locations[2]]
  s <- cov(Z)*(n-1)/n
  wi <- gtlp(s,lambda,tau)
  return(list(tau=tau,lambda=lambda,wi=wi))
}




