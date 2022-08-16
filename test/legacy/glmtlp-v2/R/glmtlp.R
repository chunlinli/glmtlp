
#' fit a GLM with l0, l1, or tlp regularization
#' 
#' Fit a generalize linear model via penalized maximum likelihood. The 
#'   regularization path is computed for the l0, lasso, or truncated lasso 
#'   penalty at a grid of values for the regularization parameter: number of 
#'   non-zeros or lambda. Can deal with all shapes of data (not large sparse 
#'   data matrices now). Fits linear, logistic, multinomial, and poisson models.
#'   
#' The sequence of models implied by \code{lambda} or \code{kappa} is fit by 
#'   coordinate descent. 
#'
#' The objective function for \code{"gaussian"} is \deqn{1/2 RSS/nobs +
#'   \lambda*penalty,} and for the other models it is \deqn{-loglik/nobs + 
#'   \lambda*penalty.} 
#'   Note also that for \code{"gaussian"}, \code{glmtlp} standardizes y to 
#'   have unit variance (using 1/(n-1) formula). The coefficients for any 
#'   predictor variables with zero variance are set to zero for all values of 
#'   lambda.
#'
#' ## Details on `family` option
#' 
#' glmtlp currently only supports built-in families.
#'
#' The built in families are specified via a character string. For all families,
#' the object produced is a regularization path for fitting the generalized 
#' linear regression paths, by maximizing the appropriate penalized
#' log-likelihood. Sometimes the sequence is truncated before \code{nlambda} 
#' values of \code{lambda} have been used, because of instabilities in the 
#' inverse link functions near a saturated fit. \code{glmtlp(..., family="binomial")} 
#' fits a traditional logistic regression model for the log-odds.
#' \code{glmtlp(..., family="multinomial")} fits a symmetric multinomial model,
#' where each class is represented by a linear model (on the log-scale). The
#' penalties take care of redundancies. A two-class \code{"multinomial"} model
#' will produce the same fit as the corresponding \code{"binomial"} model,
#' except the pair of coefficient matrices will be equal in magnitude and
#' opposite in sign, and half the \code{"binomial"} values.
#'
#' ## Details on `penalty` option
#' 
#' The built in penalties are specified by a character string. For \code{l0} 
#'   penalty, \code{kappa} sequence is used for generating the regularization 
#'   path, while for the other penalties, \code{lambda} sequence is used for 
#'   generating the regularization path.
#'
#' @param X input matrix, of dimension \code{nobs} x \code{nvars}; 
#'   each row is  an observation vector. Currently, doesn't support the sparse 
#'   matrix format (inherit from class \code{"sparseMatrix"} as in package 
#'   \code{Matrix}).
#' @param y response variable, of length nobs. Quantitative for 
#'   \code{family="gaussian"}, or \code{family="poisson"} (non-negative counts). 
#'   For \code{family="binomial"}, should be either a factor with two levels. 
#'   For \code{family="multinomial"}, (to be added). 
#' @param family a character string representing one of the built-in families. 
#'   See Details section below.
#' @param penalty a character string representing one of the built-in penalties. 
#'   "l0" represents \code{L_0} penalty, "l1" represents lasso-type penalty, and 
#'   "tlp" represents truncated lasso penalty.
#' @param nlambda The number of \code{lambda} values. Default is 100.
#' @param lambda.min.ratio Smallest value for \code{lambda}, as a fraction of 
#'   \code{lambda.max}, the (data derived) entry value (i.e. the smallest value 
#'   for which all coefficients are zero). The default depends on the sample size 
#'   \code{nobs} relative to the number of variables \code{nvars}. If \code{nobs 
#'   > nvars}, the default is \code{0.0001}, close to zero.  If \code{nobs < 
#'   nvars}, the default is \code{0.01}.
#' @param lambda A user supplied \code{lambda} sequence. Typical usage is to
#'   have the program compute its own \code{lambda} sequence based on 
#'   \code{nlambda} and \code{lambda.min.ratio}. Supplying a value of 
#'   \code{lambda} overrides this. WARNING: use with care. Avoid supplying a 
#'   single value for \code{lambda} (for predictions after CV use 
#'   \code{predict()} instead).  Supply instead a decreasing sequence of 
#'   \code{lambda} values. \code{glmtlp} relies on its warms starts for speed, 
#'   and its often faster to fit a whole path than compute a single fit. When
#'   the \code{penalty} is 'l0', it is not recommended for the users to supply 
#'   this parameter.
#' @param kappa A user supplied \code{kappa} sequence. Typical usage is to
#'   have the program compute its own \code{kappa} sequence based on \code{nvars} 
#'   and \code{nobs}. This sequence is used when penalty is 'l0'.
#' @param tau A tuning parameter used in the TLP penalty. Default is 
#'   \code{0.3 * sqrt(log(nvars)/nobs)}.
#' @param delta A tuning parameter used in the coordinate majorization descent 
#'   algorithm. See Yang, Y., & Zou, H. (2014) in the reference for detail.
#' @param tol Tolerance level for all the iterative optimization algorithms.
#' @param weights observation weights. Default is 1 for each observation
#' @param penalty.factor Separate penalty factors can be applied to each 
#'   coefficient. This is a number that multiplies \code{lambda} to allow 
#'   differential shrinkage. Can be 0 for some variables, which implies no 
#'   shrinkage, and that variable is always included in the model. Default is 1 
#'   for all variables. 
#' @param standardize Logical. Whether or not standardize the input matrix 
#'   \code{X}; default is \code{TRUE}.
#' @param dc.maxit Maximum number of iterations for the DC (Difference of 
#'   Convex Functions) programming; default is 20.
#' @param cd.maxit Maximum number of iterations for the coordinate descent 
#'   algorithm; default is 10^4.
#' @param nr.maxit Maximum number of iterations for the Newton-Raphson method; 
#'   default is 500.
#' @param ... Additional argument. These include some of the original arguments 
#'   to 'glmtlp', and each must be named if used.
#' @return An object with S3 class \code{"glmtlp"}. 
#' 
#' \item{beta}{a \code{nvars x length(kappa)} matrix of
#'   coefficients when penalty is 'l0'; a \code{nvars x length(lambda)} matrix of
#'   coefficients when penalty is 'l1' and 'tlp'.}
#' \item{call}{the call that produced this object.} 
#' \item{family}{the distribution family used in the model fitting.}
#' \item{intercept}{the intercept vector, of length(kappa) or length(lambda).}
#' \item{lambda}{The actual sequence of \code{lambda} values used. Note that 
#'   the length may be smaller than the provided \code{nlambda} due to removal 
#'   of saturated values.} 
#' \item{penalty}{the penalty type in the model fitting.}
#' \item{penalty.factor}{the penalty factor for each coefficient used in the model fitting.}
#' \item{tau}{the tuning parameter used in the model fitting, available when 
#'   penalty = 'tlp'.}
#' \item{user.lambda}{logical, whether or not the used \code{lambda} sequence is 
#'   provided by the user.}
#' \item{user.kappa}{logical, whether or not the used \code{kappa} sequence is 
#'   provided by the user, available when penalty = 'l0'.}
#' 
#' @author Chunlin Li, Yu Yang, Chong Wu
#'   \cr Maintainer: Yu Yang \email{yang6367@umn.edu}
#' 
#' @seealso \code{print}, \code{predict}, \code{coef} and \code{plot} methods,
#' and the \code{cv.glmtlp} function.
#' 
#' @references Shen, X., Pan, W., & Zhu, Y. (2012). 
#'   \emph{Likelihood-based selection and sharp parameter estimation. 
#'   Journal of the American Statistical Association, 107(497), 223-232.}
#'   \cr Yang, Y., & Zou, H. (2014). \emph{A coordinate majorization descent algorithm 
#'   for l1 penalized learning. Journal of Statistical Computation and 
#'   Simulation, 84(1), 84-95.}
#'   \cr Two R package Github: \emph{ncvreg} and \emph{glmnet}.
#'   
#' @keywords models regression
#' 
#' @examples
#'
#' # Gaussian
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' fit1 <- glmtlp(X, y, family = "gaussian", penalty = "l0")
#' fit2 <- glmtlp(X, y, family = "gaussian", penalty = "l1")
#' fit3 <- glmtlp(X, y, family = "gaussian", penalty = "tlp")
#' 
#' # Binomial
#' 
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- sample(c(0,1), 100, replace = TRUE)
#' fit <- glmtlp(X, y, family = "binomial", penalty = "l1")
#'
#' @importFrom stats model.matrix
#' @export glmtlp

glmtlp <- function(X, y, family=c("gaussian","binomial"), penalty=c("l0", "l1", "tlp"),
                   nlambda=100, lambda.min.ratio=ifelse(nobs < nvars, 1e-2, 1e-4), 
                   lambda=NULL, kappa=NULL, 
                   tau=0.3*sqrt(log(nvars)/nobs), delta=2.0, tol=1e-4, 
                   weights=NULL, penalty.factor=rep(1.0, nvars), standardize=TRUE, 
                   dc.maxit=20, cd.maxit=10000, nr.maxit=500, ...) {
  # Coersion
  this.call <- match.call()
  family <- match.arg(family)
  penalty <- match.arg(penalty)

  # check X
  xdim <- dim(X)
  if (is.null(xdim) | (xdim[2] <= 1)) stop("X should be a matrix with 2 or more columns")
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (typeof(X)=="character") stop("X must be a numeric matrix")
  if (typeof(X)=="integer") storage.mode(X) <- "double"
  
  nobs <- as.integer(xdim[1])
  nvars <- as.integer(xdim[2])

  # check y
  y <- drop(y) # Delete the dimensions of an array which has only one level.
  dimy <- dim(y)
  nrowy <- ifelse(is.null(dimy), length(y), dimy[1])
  if (nrowy != nobs) stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of X (", nobs, ")", sep = ""))
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric"),
      as.double(y))
    options(op)
  }
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to the model")
  if (family == "binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data", call.=FALSE)
  if (family == "binomial" & !identical(sort(unique(y)), 0:1)) {
    y <- as.double(y == max(y))
  } 

  # check penalty.factor, weights, and s
  penalty.factor <- as.double(penalty.factor) 
  if (length(penalty.factor) != nvars) stop(paste("the length of penalty.factor (", length(penalty.factor), ") not equal to the number of variables (", nvars, ")", sep = ""))
  if (is.null(weights)) {
      weights <- rep(1.0, nobs)
  } else if (length(weights) != nobs) {
      stop(paste("number of elements in weights (", length(weights), ") not equal to the number of rows of X (", nobs, ")", sep = ""))
  } else {
      weights <- as.double(weights)
  }
  
  ## Deprication support
  dots <- list(...)
#   if ("n.lambda" %in% names(dots)) nlambda <- dots$n.lambda

  ## standardize X
  if (standardize) {
    X <- std(X)
  }

  ## check/setup lambda and kappa
  if (is.null(lambda)) {
    nlambda <- as.integer(nlambda)
    if (nlambda < 2) stop("nlambda must be at least 2")
    if(lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
    lambda <- setup_lambda(X, y, lambda.min.ratio, nlambda)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    if (nlambda < 2) stop("the length of lambda must be at least 2")
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    if (!is.double(lambda)) lambda <- as.double(lambda)
    lambda <- sort(lambda, decreasing = TRUE)
    user.lambda <- TRUE
  }

  if (penalty == "l0") {
    if(is.null(kappa)) {
      kappa <- 1:min(nvars, as.integer(nvars/log(nobs)))
      nkappa <- length(kappa)
      user.kappa <- FALSE
    } else {
      if (any(kappa < 0)) stop("kappa should be non-negative")
      if (!is.integer(kappa)) {
        message("Coercing kappa to integers.")
        kappa <- as.integer(kappa) 
      }
      kappa <- sort(kappa, decreasing = FALSE)
      nkappa <- length(kappa)
      user.kappa <- TRUE
    }
  }
  
  ## check tau, delta, tol, and maxiters: may add more on tau, delta, and tol checking
  tau <- as.double(tau)
  delta <- as.double(delta)
  tol <- as.double(tol)
  dc.maxit <- as.integer(dc.maxit)
  cd.maxit <- as.integer(cd.maxit)
  nr.maxit <- as.integer(nr.maxit)
  
  ## fit
  if (family == "gaussian") {
    fit <- switch(penalty,
                  "l0" = gaussian_l0_reg(X, y, kappa, lambda, tau, weights, penalty.factor, delta, tol, dc.maxit, cd.maxit), 
                  "l1" = gaussian_l1_reg(X, y, lambda, weights, penalty.factor, delta, tol, cd.maxit), 
                  "tlp" = gaussian_tlp_reg(X, y, lambda, tau, weights, penalty.factor, delta, tol, dc.maxit, cd.maxit)
                  ) 
  } else if (family == "binomial") {
    fit <- switch(penalty,
                  "l0" = logistic_l0_reg(), 
                  "l1" = logistic_l1_reg(X, y, lambda, penalty.factor, delta, tol, nr.maxit, cd.maxit), 
                  "tlp" = logistic_tlp_reg()
                  )
  }

  ## Unstandardize
  # beta <- matrix(0, nrow = nvars + 1, ncol = ifelse(penalty == "l0", nkappa, nlambda))
  # b0 <- fit$b0
  # b <- fit$b
  # bb <- b / attr(X, "scaled:scale")
  # beta[2:(nvars+1), ] <- bb
  # beta[1, ] <- b0 - crossprod(attr(X, "scaled:center"), bb)
  
  intercept <- fit$b0
  beta <- fit$b
  
  if (standardize) {
    beta <- beta / attr(X, "scaled:scale")
    intercept <- intercept - crossprod(attr(X, "scaled:center"), beta)
  } 
  

  ## Names
  varnames <- colnames(X)
  if (is.null(varnames)) varnames <- paste("V", seq(nvars), sep = "") 
  rownames(beta) <- varnames
  if (penalty == "l0") {
    colnames(beta) <- paste(kappa)
  } else {
    colnames(beta) <- lambda_names(lambda)
  }
  ## Output
  out <- structure(list(beta = beta, 
                        call = this.call,
                        family = family,
                        intercept = intercept, 
                        penalty = penalty, 
                        penalty.factor = penalty.factor, 
                        weights = weights
                        ), 
                   class = "glmtlp")
  if (penalty == "l0") {
    out$lambda <- lambda
    out$user.lambda <- user.lambda
    out$kappa <- kappa
    out$user.kappa <- user.kappa
  } else if (penalty == "l1") {
    out$lambda <- lambda
    out$user.lambda <- user.lambda
  } else {
    out$lambda <- lambda
    out$user.lambda <- user.lambda
    out$tau <- tau
  }
  out
}

