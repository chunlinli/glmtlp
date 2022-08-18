#' Solve Truncated L1 Constrained Problem
#'
#' Using Difference Convex (DC) method to solve minimization with Truncated L1 Constraint.
#'
#' @param Y Response vector of length n.
#' @param X Design matrix with n rows and p columns.
#' @param K Constraint parameter in TLC.
#' @param tau Parameter tau in TLC.
#' @param tlc_weight Length p vector of weights corresponding to p variables,
#' each element is either 1 (with constraint) or 0 (without constraint).
#' Default are all 1's.
#' @param maxit_tlc Maximum number of DC iteration.
#' @param tol_tlc Convergence tolerance for TLC.
#'
#' @return beta, length p vector of estimated coefficients.
#' @export
#'
#' @examples
TLC <- function(Y,X,K,
                tau = 1e-5,
                tlc_weight = rep(1,ncol(X)),
                maxit_tlc = 100,tol_tlc=1e-5)
{
  n = length(Y)
  p = ncol(X)

  Y = scale(Y,scale = F)
  X = scale(X,scale = F)

  beta = rep(0,p)
  loss_old = sum((Y - X%*%beta)^2)
  loss_new = loss_old

  ite_ind = 1
  while(ite_ind <= maxit_tlc)
  {
    ite_ind = ite_ind+1
    weight_lasso = (tlc_weight / tau) * (abs(beta) <= tau)
    t = K - sum(tlc_weight * (abs(beta) > tau))
    if(t<=0)
    {
      zero_ind = which(weight_lasso==0)
      if(length(zero_ind) == 0)
      {
        lm1 = lm(Y ~ 0)
      } else {
        lm1 = lm(Y ~ 0 + X[,zero_ind])
      }
      beta1 = lm1$coefficients
      names(beta1) = NULL
      beta = rep(0,p)
      beta[zero_ind] = beta1
    } else {
      weight_lasso = weight_lasso*tau
      t = t*tau
      zero_ind = which(weight_lasso==0)
      data1 = as.data.frame(cbind(Y,X))
      fit_formula = as.formula(paste(colnames(data1)[1], "~."))

      if(length(zero_ind > 0))
      {
        sweep_formula = as.formula(paste("~1+",
                                         paste(colnames(data1)[zero_ind+1],
                                               collapse = "+"),sep=""))

      } else {
        sweep_formula = as.formula(paste("~1"))
      }
      lasso2fit = lasso2::l1ce(fit_formula,data1,bound=t,absolute.t = T,
                               sweep.out = sweep_formula,standardize = F)
      beta = lasso2fit$coefficients[-1]
      names(beta) = NULL
    }
    loss_new = sum((Y - X%*%beta)^2)
    if(abs(loss_new - loss_old) < tol_tlc)
    {
      break()
    } else{
      loss_old = loss_new
    }
  }
  nonzero_ind = which(beta!=0)
  if(length(nonzero_ind) == 0)
  {
    lm1 = lm(Y ~ 0)
  } else{
    lm1 = lm(Y ~ 0 + X[,nonzero_ind])
  }
  beta1 = lm1$coefficients
  names(beta1) = NULL
  beta = rep(0,p)
  beta[nonzero_ind] = beta1

  return(beta)

}
