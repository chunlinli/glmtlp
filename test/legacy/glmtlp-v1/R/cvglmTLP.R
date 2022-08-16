cv.glmTLP <- function(x, y, family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),nfolds = 10, weights, offset=NULL, lambda, tau = 0.3, nlambda=100, penalty.factor = rep(1, nvars),  lambda.min.ratio=ifelse(nobs<nvars,1e-3,1e-4),standardize=TRUE,intercept=TRUE,dfmax=nvars+1,pmax=min(dfmax*2+20,nvars),lower.limits=-Inf,upper.limits=Inf,standardize.response=FALSE, maxIter=100, Tol=1e-4){
    
    eps0=1e-3
    
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
    
    
    folds.i <- sample(rep(1:nfolds, length.out = nobs))
    
    cv.tmp <- matrix(NA, nrow = nfolds, ncol = length(lambda))
    for (k in 1:nfolds) {
        test.i <- which(folds.i == k)
        x.train <- x[-test.i,]
        y.train <- y[-test.i]
        
        weights.tmp <- weights[-test.i]
        fitobj <- glmTLP(x.train, y.train, family=family, weights = weights.tmp, offset=offset, lambda=lambda, penalty.factor = penalty.factor, tau = tau, standardize = standardize, intercept = intercept, dfmax=dfmax, pmax=pmax, lower.limits = lower.limits, upper.limits = upper.limits, standardize.response = standardize.response, maxIter = maxIter, Tol = Tol)

        x.test <- x[test.i,]
        y.test <- y[test.i]
        
        preds <- predict(fitobj,x.test,s = lambda,"response")
        
        cv.tmp[k, ] <- colSums((y.test-preds)^2)
    }
    cvm <- colMeans(cv.tmp)
    cvsd <- apply(cv.tmp, 2, sd)
    cvup <- apply(cv.tmp, 2, max)
    cvlo <- apply(cv.tmp, 2, min)
    
    glmnet.fit <- glmTLP(x, y, family=family, weights = weights, offset=offset, lambda=lambda, penalty.factor = penalty.factor, tau = tau, standardize = standardize, intercept = intercept,  dfmax=dfmax, pmax=pmax, lower.limits = lower.limits, upper.limits = upper.limits, standardize.response = standardize.response, maxIter = maxIter, Tol = Tol)
    
    lambda.min <- lambda[cvm ==min(cvm)]
    lambda.1se <- lambda[cvm <= (min(cvm) + sd(cvm)) ]
    lambda.1se <- min(lambda.1se)
    
    out <- list(lambda =lambda, tau = tau, cvm = cvm, cvsd = cvsd, cvup = cvup, cvlo = cvlo,  name = "Mean-Squared Error", glmnet.fit = glmnet.fit,  lambda.min =  lambda.min, lambda.1se = lambda.1se)
    class(out) = "cv.glmnet"
    out
}