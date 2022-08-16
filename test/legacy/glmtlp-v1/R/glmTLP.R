glmTLP<-function(x, y, family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),weights, offset=NULL, lambda, tau = 0.3, nlambda=100, penalty.factor = rep(1, nvars), lambda.min.ratio=ifelse(nobs<nvars,1e-3,1e-4),standardize=TRUE,intercept=TRUE,dfmax=nvars+1,pmax=min(dfmax*2+20,nvars),lower.limits=-Inf,upper.limits=Inf,standardize.response=FALSE, maxIter=100, Tol=1e-4){
    
    this.call=match.call()
    eps0=1e-4

    np = dim(x)
    if(is.null(np)|(np[2]<=1)) stop("x should be a matrix with 2 or more columns")
    nobs <- as.integer(np[1])
    nvars <- as.integer(np[2])
    
    ne <- as.integer(dfmax)
    nx <- as.integer(pmax)
    dfmax <- ne
    pmax <- nx
    
    family <- match.arg(family)

    if(missing(weights)) weights <- rep(1,nobs)
    else if(length(weights)!=nobs) stop(paste("number of elements in weights (",length(weights),") not equal to the number of rows of x (",nobs,")",sep=""))
    
    nlam <- as.integer(nlambda)
 
    if(missing(lambda)) {
        if(lambda.min.ratio>=1) stop("lambda.min.ratio should be less than 1")
        lambda.max <- max(abs(t(x) %*% y)) / nobs
        lambda <- exp(seq(log(lambda.max),log(lambda.min.ratio*lambda.max),len=nlam))
    } else {
        if(any(lambda<0))stop("lambdas should be non-negative")
        lambda <- sort(lambda,decreasing= TRUE)
    }
    
    
    res0 <- glmnet(x, y, family=family, weights = weights, offset=offset, lambda=lambda, penalty.factor = penalty.factor, standardize = standardize, intercept = intercept,  dfmax=dfmax, pmax=pmax, lower.limits = lower.limits, upper.limits = upper.limits, standardize.response = standardize.response)
    
    bs0 <- bs1 <- bs1.1 <- res0$beta
    a0 <- a0.1 <- res0$a0
    
    df.tmp <- rep(0,length(lambda))
    dev.ratio.tmp <- rep(0,length(lambda))
    
    npasses <- 0
    for(k in 1:length(lambda)){
        b0 <- bs0[,k]
        lambda.wt <- (abs(b0)<=tau)
        if (length(lambda.wt[abs(lambda.wt)<eps0])>nobs/2)
            lambda.wt[abs(lambda.wt)<eps0] <- eps0

        lambda.wt  <- lambda.wt * penalty.factor
        lambda.tmp <- lambda[k] * (sum(lambda.wt) / sum(penalty.factor))
        
        res <- glmnet(x, y, family=family, weights = weights, offset=offset, lambda=lambda.tmp, penalty.factor=lambda.wt, standardize = standardize, intercept = intercept,lower.limits = lower.limits, upper.limits = upper.limits, standardize.response = standardize.response)
        b1 <- bs1.1[,k] <- res$beta[,1]
        a0.1[k] <- res$a0[1]
        ##iterate:
        iter <- 1
        while(iter <= maxIter && any(abs(b1-b0)>Tol)){
            iter <- 1+iter
            b0 <- b1
            lambda.wt <- (abs(b0)<=tau)
            if (length(lambda.wt[abs(lambda.wt)<eps0])>nobs/2)
                lambda.wt[abs(lambda.wt)<eps0] <- eps0
                
            lambda.wt  <- lambda.wt * penalty.factor
            lambda.tmp <- lambda[k] * (sum(lambda.wt) / sum(penalty.factor))
            
            res <- glmnet(x, y, family=family, weights = weights, offset=offset, lambda=lambda.tmp, penalty.factor=lambda.wt, standardize = standardize, intercept = intercept,lower.limits = lower.limits, upper.limits = upper.limits, standardize.response = standardize.response)

            b1 <- res$beta[,1]
        }

        bs1[,k] <- res$beta[,1]
        a0[k] <- res$a0[1]
        df.tmp[k] <- res$df
        dev.ratio.tmp[k] <-res$dev.ratio
        npasses <- npasses + res$npasses
    }
    res0$call <- this.call
    res0$beta <- bs1
    res0$a0 <- a0
    res0$npasses <- npasses
    res0$df <- df.tmp
    res0$dev.ratio <- dev.ratio.tmp
    res0$iter <- iter
    return(res0)
}

