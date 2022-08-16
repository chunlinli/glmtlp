\name{glmTLP}
\alias{glmTLP}
\title{Fit a GLM with truncated lasso regularization}
\description{Fit a generalized linear model via penalized maximum likelihood. The regularization path is computed for the truncated lasso at a grid of values for the regularization parameter lambda. Can deal with all shapes of data, including very large sparse data matrices. Fit linear, logistic and multinomial, poisson, and Cox regression models.
}
\usage{
glmTLP(x, y, family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),
    weights, offset=NULL, lambda, tau = 0.3, nlambda=100,
    penalty.factor = rep(1, nvars), lambda.min.ratio=ifelse(nobs<nvars,1e-3,1e-4),
    standardize=TRUE,intercept=TRUE,dfmax=nvars+1,pmax=min(dfmax*2+20,nvars),
    lower.limits=-Inf,upper.limits=Inf,
    standardize.response=FALSE, maxIter=100, Tol=1e-4)
}

\arguments{
  \item{x}{input matrix, of dimension nobs x nvars; each row is an
  observation vector. Can be in sparse matrix format (inherit from class \code{"sparseMatrix"} as in package \code{Matrix}; not yet available for \code{family="cox"})}
  \item{y}{response variable. Quantitative for \code{family="gaussian"},
  or \code{family="poisson"} (non-negative counts). For
  \code{family="binomial"} should be either a factor with two levels, or
  a two-column matrix of counts or proportions (the second column is
  treated as the target class; for a factor, the last level in
  alphabetical order is the target class). For
  \code{family="multinomial"}, can be a \code{nc>=2} level factor, or a
  matrix with \code{nc} columns of counts or proportions.
  For either \code{"binomial"} or \code{"multinomial"}, if \code{y} is
  presented as a vector, it will be coerced into a factor. For
  \code{family="cox"}, \code{y} should be a two-column matrix with
  columns named 'time' and 'status'. The latter is a binary variable,
  with '1' indicating death, and '0' indicating right censored. The
  function \code{Surv()} in package \pkg{survival} produces such a
  matrix. For  \code{family="mgaussian"}, \code{y} is a matrix of quantitative responses.}
  \item{family}{Response type (see above)}
  \item{weights}{observation weights. Can be total counts if responses are proportion matrices. Default is 1 for each observation}
  \item{offset}{A vector of length \code{nobs} that is included in the linear predictor (a \code{nobs x nc} matrix for the \code{"multinomial"} family). Useful for the \code{"poisson"} family (e.g. log of exposure time), or for refining a model by starting at a current fit. Default is \code{NULL}. If supplied, then values must also be supplied to the \code{predict} function.}
  \item{tau}{Write something about \code{tau}}
  \item{nlambda}{The number of \code{lambda} values - default is 100.}
    \item{penalty.factor}{Separate penalty factors can be applied to each
    coefficient. This is a number that multiplies \code{lambda} to allow
    differential shrinkage. Can be 0 for some variables, which implies
    no shrinkage, and that variable is always included in the
    model. Default is 1 for all variables (and implicitly infinity for
    variables listed in \code{exclude}). Note: the penalty factors are
  internally rescaled to sum to nvars, and the lambda sequence will
  reflect this change.}
  \item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of
    \code{lambda.max}, the (data derived) entry value (i.e. the smallest
  value for which all coefficients are zero). The default depends on the
  sample size \code{nobs} relative to the number of variables
  \code{nvars}. If \code{nobs > nvars}, the default is \code{0.0001},
  close to zero.  If \code{nobs < nvars}, the default is \code{0.01}.
  A very small value of
  \code{lambda.min.ratio} will lead to a saturated fit in the \code{nobs <
  nvars} case. This is undefined for
  \code{"binomial"} and \code{"multinomial"} models, and \code{glmnet}
  will exit gracefully when the percentage deviance explained is almost
  1.} 
  \item{lambda}{A user supplied \code{lambda} sequence. Typical usage
    is to have the 
    program compute its own \code{lambda} sequence based on
    \code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
    \code{lambda} overrides this. WARNING: use with care. Do not supply
  a single value for \code{lambda} (for predictions after CV use \code{predict()}
  instead).  Supply instead
    a decreasing sequence of \code{lambda} values. \code{glmnet} relies
  on its warms starts for speed, and its often faster to fit a whole
  path than compute a single fit.}
  \item{standardize}{Logical flag for x variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is \code{standardize=TRUE}.
  If variables are in the same units already, you might not wish to
  standardize. See details below for y standardization with \code{family="gaussian"}.}
 \item{intercept}{Should intercept(s) be fitted (default=TRUE) or set to
  zero (FALSE)}
  \item{dfmax}{Limit the maximum number of variables in the
    model. Useful for very large \code{nvars}, if a partial path is desired.}
  \item{pmax}{Limit the maximum number of variables ever to be nonzero}
  \item{lower.limits}{Vector of lower limits for each coefficient;
  default \code{-Inf}. Each
  of these must be non-positive. Can be presented as a single value
  (which will then be replicated), else a vector of length \code{nvars}}
 \item{upper.limits}{Vector of upper limits for each coefficient;
  default \code{Inf}. See \code{lower.limits}}
\item{standardize.response}{This is for the \code{family="mgaussian"}
  family, and allows the user to standardize the response variables}
  \item{maxIter}{Maximum iteration for TLP.}
  \item{Tol}{Tolerance.}
}
\details{ Write something about the details.
  }
\value{
An object that inherits from \code{glmnet}.
  \item{call}{the call that produced this object}
  \item{a0}{Intercept sequence of length \code{length(lambda)}}
  \item{beta}{For \code{"elnet"}, \code{"lognet"}, \code{"fishnet"} and \code{"coxnet"} models, a \code{nvars x
      length(lambda)} matrix of coefficients, stored in sparse column
    format (\code{"CsparseMatrix"}). For \code{"multnet"} and \code{"mgaussian"}, a list of \code{nc} such
    matrices, one for each class.}
  \item{lambda}{The actual sequence of \code{lambda} values used. }
  \item{dev.ratio}{The fraction of (null) deviance explained (for \code{"elnet"}, this
      is the R-square). The deviance calculations incorporate weights if
  present in the model. The deviance is defined to be 2*(loglike_sat -
  loglike), where loglike_sat is the log-likelihood for the saturated
  model (a model with a free parameter per observation). Hence dev.ratio=1-dev/nulldev.}
    \item{nulldev}{Null deviance (per observation). This is defined to
  be  2*(loglike_sat -loglike(Null)); The NULL model refers to the
  intercept model, except for the Cox, where it is the 0 model.}
  \item{df}{The number of nonzero coefficients for each value of
    \code{lambda}. }
  \item{dim}{dimension of coefficient matrix (ices)}
  \item{nobs}{number of observations}
  \item{npasses}{total passes over the data summed over all lambda
    values}
  \item{offset}{a logical variable indicating whether an offset was included in the model}
  \item{jerr}{error flag, for warnings and errors (largely for internal debugging).}
}


\examples{
	data("QuickStartExample")
	fit = glmTLP(x,y)
}

\references{ Xiaotong Shen , Wei Pan and Yunzhang Zhu (2012)
  \emph{Likelihood-Based Selection and Sharp Parameter Estimation}, 
  \emph{Journal of the American Statistical Association, 107:497, 223-232}}
\author{Chong Wu, Wei Pan\cr 
Maintainer: Chong Wu <wuxx0845@umn.edu>}

\keyword{models}
\keyword{regression}

 
