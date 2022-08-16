\name{cv.glmTLP}
\alias{cv.glmTLP}
\title{Cross-validation for glmTLP}
\description{Does k-fold cross-validation for glmTLP, produces a plot, and returns a value for \code{lambda} with pre-specified \code{tau}.}
\usage{
cv.glmTLP(x, y, family=c("gaussian","binomial","poisson","multinomial","cox","mgaussian"),
nfolds = 10, weights, offset=NULL, lambda, tau = 0.3, 
nlambda=100, penalty.factor = rep(1, nvars), 
 lambda.min.ratio=ifelse(nobs<nvars,1e-3,1e-4),
standardize=TRUE,intercept=TRUE,dfmax=nvars+1,
pmax=min(dfmax*2+20,nvars), lower.limits=-Inf,upper.limits=Inf,
standardize.response=FALSE, maxIter=100, Tol=1e-4)
}
\arguments{
  \item{x}{\code{x} matrix as in \code{glmnet}.}
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
    \item{nfolds}{number of folds - default is 10. Although \code{nfolds}
    can be as large as the sample size (leave-one-out CV), it is not
    recommended for large datasets. Smallest value allowable is \code{nfolds=3}}
  \item{weights}{Observation weights; defaults to 1 per observation}
  \item{offset}{Offset vector (matrix) as in \code{glmnet}}
  \item{lambda}{Optional user-supplied lambda sequence; default is
    \code{NULL}, and \code{glmTLP} chooses its own sequence}
  \item{tau}{Tuning parameter.}
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
\details{The function runs \code{glmTLP} \code{nfolds}+1 times; the
  first to get the \code{lambda} sequence, and then the remainder to
  compute the fit with each of the folds omitted. The error is
  accumulated, and the average error and standard deviation over the
  folds is computed. 
  Note that \code{cv.glmnet} does NOT search for
  values for \code{tau}. A specific value should be supplied, else
  \code{tau= 0.3} is assumed by default.
}
\value{an object of class \code{"cv.glmnet"} is returned, which is a
  list with the ingredients of the cross-validation fit. Although the implementation is different, we try to mimic returning as \code{"cv.glment"} in a popular package \code{glmnet} such that users can use truncated lasso as using elastic net.
  \item{lambda}{the values of \code{lambda} used in the fits.}
  \item{cvm}{The mean cross-validated error - a vector of length
    \code{length(lambda)}.}
  \item{cvsd}{estimate of standard error of \code{cvm}.}
  \item{cvup}{upper curve = \code{cvm+cvsd}.}
  \item{cvlo}{lower curve = \code{cvm-cvsd}.}
  \item{nzero}{number of non-zero coefficients at each \code{lambda}.}
  \item{name}{a text string indicating type of measure (for plotting
    purposes).}
 \item{glmnet.fit}{a fitted glmnet object for the full data.}
  \item{lambda.min}{value of \code{lambda} that gives minimum
    \code{cvm}.}
  \item{lambda.1se}{largest value of \code{lambda} such that error is
    within 1 standard error of the minimum.}
  \item{fit.preval}{if \code{keep=TRUE}, this is the array of
    prevalidated fits. Some entries can be \code{NA}, if that and
    subsequent values of \code{lambda} are not reached for that fold}
  \item{foldid}{if \code{keep=TRUE}, the fold assignments used}
}

\examples{
	data("QuickStartExample")
	fit = cv.glmTLP(x,y,tau = 1)
}

\references{ Xiaotong Shen , Wei Pan and Yunzhang Zhu (2012)
  \emph{Likelihood-Based Selection and Sharp Parameter Estimation}, 
  \emph{Journal of the American Statistical Association, 107:497, 223-232}
}
\author{Chong Wu\cr 
Maintainer: Chong Wu <wuxx0845@umn.edu>}

\keyword{models}
\keyword{regression}

