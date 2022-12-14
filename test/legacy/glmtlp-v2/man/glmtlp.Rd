% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmtlp.R, R/glmtlp_package.R
\docType{package}
\name{glmtlp}
\alias{glmtlp}
\title{fit a GLM with l0, l1, or tlp regularization}
\usage{
glmtlp(
  X,
  y,
  family = c("gaussian", "binomial"),
  penalty = c("l0", "l1", "tlp"),
  nlambda = 100,
  lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL,
  kappa = NULL,
  tau = 0.3 * sqrt(log(nvars)/nobs),
  delta = 2,
  tol = 1e-04,
  weights = NULL,
  penalty.factor = rep(1, nvars),
  standardize = TRUE,
  dc.maxit = 20,
  cd.maxit = 10000,
  nr.maxit = 500,
  ...
)
}
\arguments{
\item{X}{input matrix, of dimension \code{nobs} x \code{nvars}; 
each row is  an observation vector. Currently, doesn't support the sparse 
matrix format (inherit from class \code{"sparseMatrix"} as in package 
\code{Matrix}).}

\item{y}{response variable, of length nobs. Quantitative for 
\code{family="gaussian"}, or \code{family="poisson"} (non-negative counts). 
For \code{family="binomial"}, should be either a factor with two levels. 
For \code{family="multinomial"}, (to be added).}

\item{family}{a character string representing one of the built-in families. 
See Details section below.}

\item{penalty}{a character string representing one of the built-in penalties. 
"l0" represents \code{L_0} penalty, "l1" represents lasso-type penalty, and 
"tlp" represents truncated lasso penalty.}

\item{nlambda}{The number of \code{lambda} values. Default is 100.}

\item{lambda.min.ratio}{Smallest value for \code{lambda}, as a fraction of 
\code{lambda.max}, the (data derived) entry value (i.e. the smallest value 
for which all coefficients are zero). The default depends on the sample size 
\code{nobs} relative to the number of variables \code{nvars}. If \code{nobs 
> nvars}, the default is \code{0.0001}, close to zero.  If \code{nobs < 
nvars}, the default is \code{0.01}.}

\item{lambda}{A user supplied \code{lambda} sequence. Typical usage is to
have the program compute its own \code{lambda} sequence based on 
\code{nlambda} and \code{lambda.min.ratio}. Supplying a value of 
\code{lambda} overrides this. WARNING: use with care. Avoid supplying a 
single value for \code{lambda} (for predictions after CV use 
\code{predict()} instead).  Supply instead a decreasing sequence of 
\code{lambda} values. \code{glmtlp} relies on its warms starts for speed, 
and its often faster to fit a whole path than compute a single fit. When
the \code{penalty} is 'l0', it is not recommended for the users to supply 
this parameter.}

\item{kappa}{A user supplied \code{kappa} sequence. Typical usage is to
have the program compute its own \code{kappa} sequence based on \code{nvars} 
and \code{nobs}. This sequence is used when penalty is 'l0'.}

\item{tau}{A tuning parameter used in the TLP penalty. Default is 
\code{0.3 * sqrt(log(nvars)/nobs)}.}

\item{delta}{A tuning parameter used in the coordinate majorization descent 
algorithm. See Yang, Y., & Zou, H. (2014) in the reference for detail.}

\item{tol}{Tolerance level for all the iterative optimization algorithms.}

\item{weights}{observation weights. Default is 1 for each observation}

\item{penalty.factor}{Separate penalty factors can be applied to each 
coefficient. This is a number that multiplies \code{lambda} to allow 
differential shrinkage. Can be 0 for some variables, which implies no 
shrinkage, and that variable is always included in the model. Default is 1 
for all variables.}

\item{standardize}{Logical. Whether or not standardize the input matrix 
\code{X}; default is \code{TRUE}.}

\item{dc.maxit}{Maximum number of iterations for the DC (Difference of 
Convex Functions) programming; default is 20.}

\item{cd.maxit}{Maximum number of iterations for the coordinate descent 
algorithm; default is 10^4.}

\item{nr.maxit}{Maximum number of iterations for the Newton-Raphson method; 
default is 500.}

\item{...}{Additional argument. These include some of the original arguments 
to 'glmtlp', and each must be named if used.}
}
\value{
An object with S3 class \code{"glmtlp"}. 

\item{beta}{a \code{nvars x length(kappa)} matrix of
  coefficients when penalty is 'l0'; a \code{nvars x length(lambda)} matrix of
  coefficients when penalty is 'l1' and 'tlp'.}
\item{call}{the call that produced this object.} 
\item{family}{the distribution family used in the model fitting.}
\item{intercept}{the intercept vector, of length(kappa) or length(lambda).}
\item{lambda}{The actual sequence of \code{lambda} values used. Note that 
  the length may be smaller than the provided \code{nlambda} due to removal 
  of saturated values.} 
\item{penalty}{the penalty type in the model fitting.}
\item{penalty.factor}{the penalty factor for each coefficient used in the model fitting.}
\item{tau}{the tuning parameter used in the model fitting, available when 
  penalty = 'tlp'.}
\item{user.lambda}{logical, whether or not the used \code{lambda} sequence is 
  provided by the user.}
\item{user.kappa}{logical, whether or not the used \code{kappa} sequence is 
  provided by the user, available when penalty = 'l0'.}
}
\description{
Fit a generalize linear model via penalized maximum likelihood. The 
  regularization path is computed for the l0, lasso, or truncated lasso 
  penalty at a grid of values for the regularization parameter: number of 
  non-zeros or lambda. Can deal with all shapes of data (not large sparse 
  data matrices now). Fits linear, logistic, multinomial, and poisson models.
  
The sequence of models implied by \code{lambda} or \code{kappa} is fit by 
  coordinate descent.

The package provides 3 penalties: l0, l1, and tlp and 3 distribution families: 
  gaussian, binomial, and poisson.
}
\details{
The objective function for \code{"gaussian"} is \deqn{1/2 RSS/nobs +
  \lambda*penalty,} and for the other models it is \deqn{-loglik/nobs + 
  \lambda*penalty.} 
  Note also that for \code{"gaussian"}, \code{glmtlp} standardizes y to 
  have unit variance (using 1/(n-1) formula). The coefficients for any 
  predictor variables with zero variance are set to zero for all values of 
  lambda.

## Details on `family` option

glmtlp currently only supports built-in families.

The built in families are specified via a character string. For all families,
the object produced is a regularization path for fitting the generalized 
linear regression paths, by maximizing the appropriate penalized
log-likelihood. Sometimes the sequence is truncated before \code{nlambda} 
values of \code{lambda} have been used, because of instabilities in the 
inverse link functions near a saturated fit. \code{glmtlp(..., family="binomial")} 
fits a traditional logistic regression model for the log-odds.
\code{glmtlp(..., family="multinomial")} fits a symmetric multinomial model,
where each class is represented by a linear model (on the log-scale). The
penalties take care of redundancies. A two-class \code{"multinomial"} model
will produce the same fit as the corresponding \code{"binomial"} model,
except the pair of coefficient matrices will be equal in magnitude and
opposite in sign, and half the \code{"binomial"} values.

## Details on `penalty` option

The built in penalties are specified by a character string. For \code{l0} 
  penalty, \code{kappa} sequence is used for generating the regularization 
  path, while for the other penalties, \code{lambda} sequence is used for 
  generating the regularization path.
}
\section{glmtlp functions}{

`glmtlp()`, `cv.glmtlp()`
}

\examples{

# Gaussian
X <- matrix(rnorm(100 * 20), 100, 20)
y <- rnorm(100)
fit1 <- glmtlp(X, y, family = "gaussian", penalty = "l0")
fit2 <- glmtlp(X, y, family = "gaussian", penalty = "l1")
fit3 <- glmtlp(X, y, family = "gaussian", penalty = "tlp")

# Binomial

X <- matrix(rnorm(100 * 20), 100, 20)
y <- sample(c(0,1), 100, replace = TRUE)
fit <- glmtlp(X, y, family = "binomial", penalty = "l1")

}
\references{
Shen, X., Pan, W., & Zhu, Y. (2012). 
  \emph{Likelihood-based selection and sharp parameter estimation. 
  Journal of the American Statistical Association, 107(497), 223-232.}
  \cr Yang, Y., & Zou, H. (2014). \emph{A coordinate majorization descent algorithm 
  for l1 penalized learning. Journal of Statistical Computation and 
  Simulation, 84(1), 84-95.}
  \cr Two R package Github: \emph{ncvreg} and \emph{glmnet}.
}
\seealso{
\code{print}, \code{predict}, \code{coef} and \code{plot} methods,
and the \code{cv.glmtlp} function.
}
\author{
Chunlin Li, Yu Yang, Chong Wu
  \cr Maintainer: Yu Yang \email{yang6367@umn.edu}
}
\keyword{models}
\keyword{regression}
