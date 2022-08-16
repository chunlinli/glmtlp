source("/Users/lichen/Dropbox/Inference for a causal interventional model with confounders/code/tlpreg.r",chdir = TRUE)

library(mnormt)

cv.MB_Union <- function(R,tau.list,gamma.list,n.fold){
  p <- ncol(R)
  taus <- numeric(p)
  gammas <- numeric(p)
  for (i in 1:p){
    res <- cv.tlpreg0.aic(R[,i],R[,-i,drop=FALSE],tau.list,gamma.list,n.fold)
    taus[i] <- res$tau
    gammas[i] <- res$gamma
  }
  S = MB_Union(R,taus,gammas)
  return(list(S=S,taus=taus,gammas=gammas))
}

MB_Union <- function(R,taus,gammas){
  p <- ncol(R)
  q <- p-1
  V <- matrix(0,p,p)
  diag(V) <- 1
  for (i in 1:p){
    y <- R[,i]
    X <- R[,-i,drop=FALSE]
    b <- tlpreg0(y,X=X,tau=taus[i],gamma=gammas[i])
    ind <- b!=0
    b.ols <- matrix(0,q,1)
    if (sum(ind)==0){
      b.ols[,1] <- b
    }else{
      aic.memo <- numeric(sum(ind))
      ord <- order(abs(b),decreasing = TRUE)
      for (l in 1:sum(ind)){
        ind2 <- which(abs(b)>=abs(b[ord[l]]))
        mod <- lm(y~X[,ind2,drop=FALSE])
        aic.memo[l] <- AIC(mod)
      }
      l <- which.min(aic.memo)
      ind <- which(abs(b)>=abs(b[ord[l]]))
      mod <- lm(y~X[,ind,drop=FALSE])
      b.ols[ind,1] <- as.numeric(mod$coefficients)[-1]
    }
    V[-i,i] <- b.ols[,1]
  }
  
  S <- 1*(V*t(V)!=0)
  return(S)
}

precision_refit <- function(Sigma,S,max.it,tol){
  
  p <- nrow(Sigma)
  n.nonzero <- sum(as.vector(S)==1)
  n.nonzero <- p+(n.nonzero-p)/2
  V1 <- vector("list",length = n.nonzero)
  V2 <- vector("list",length = n.nonzero)
  errorflag <- 0
  
  nonzero.flag <- 1
  for (i in 1:p){
    for (j in i:p){
      if (S[i,j] == 1){
        V <- matrix(0,p,p)
        V[i,j] <- 1
        V[j,i] <- 1
        V1[[nonzero.flag]] <- V
        V2[[nonzero.flag]] <- sum(diag(V%*%Sigma))
        nonzero.flag <- nonzero.flag + 1
      }
    }
  }
  
  omega.init <- diag(p)
  omega.prev <- omega.curr <- omega.init
  T.prev <- T.curr <- -log(det(omega.init)) + sum(diag(omega.init%*%Sigma))
  tol.iter <- 1
  iter <- 1
  while(tol.iter>tol&iter<max.it){
    for(i in 1:n.nonzero){
      t <- step.finding(omega.curr,V1[[i]],V2[[i]])
      omega.curr <- omega.curr+t*V1[[i]]
    }
    T.curr <- -log(det(omega.curr)) + sum(diag(omega.curr%*%Sigma))
    tol.iter <- T.prev - T.curr
    if(tol.iter<0){
      omega.curr <- omega.prev
      errorflag <- 1
      break;
    }
    omega.prev <- omega.curr
    T.prev <- T.curr
    iter <- iter+1
  }
  
  if(iter==max.it){
    warning("Coordinate Descent Algorithm for precision matrix refit does not converge!")
  }
  if(errorflag==1){
    warning("Computational issue: Ascent occurred when minimizing loss!")
  }
  return(omega.curr)
}

step.finding <- function(omega,V,c){
  U <- pd.solve(omega)
  Z <- U%*%V
  e2 <- eigen(Z)
  lambdas <- e2$values
  f0 <- f(0,lambdas,c)
  if(f0==0){
    t <- 0
  }else if(f0<0){
    if(min(lambdas)<0){
      cei <- -1/min(lambdas)
      up <- (cei)/2
      fup <- f(up,lambdas,c)
      while(fup<0){
        up <- (up+cei)/2
        fup <- f(up,lambdas,c)
      }
      t <- uniroot(f,lambdas=lambdas,c=c,interval = c(0,up))$root
    }else{
      up <- 1
      fup <- f(up,lambdas,c)
      while(fup<0){
        up <- up*2
        fup <- f(up,lambdas,c)
      }
      t <- uniroot(f,lambdas=lambdas,c=c,interval = c(0,up))$root
    }
  }else{
    if(max(lambdas)>0){
      flo <- -1/max(lambdas)
      low <- flo/2
      flow <- f(low,lambdas,c)
      while(flow>0){
        low <- (flo+low)/2
        flow <- f(low,lambdas,c)
      }
      t <- uniroot(f,lambdas=lambdas,c=c,interval = c(low,0))$root
    }else{
      low <- -1
      flow <- f(low,lambdas,c)
      while(flow>0){
        low <- 2*low
        flow <- f(low,lambdas,c)
      }
      t <- uniroot(f,lambdas=lambdas,c=c,interval = c(low,0))$root
    }
  }
  return(t)
}

f <- function(t,lambdas,c){
  l <- length(lambdas)
  result <- -sum(lambdas/(1+t*lambdas))+c
  return(result)
}


