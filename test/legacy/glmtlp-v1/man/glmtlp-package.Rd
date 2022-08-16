\name{glmtlp-package}
\alias{glmtlp-package}
\docType{package}
\title{
Truncated lasso model paths for some generalized linear models
}
\description{
  This package fits truncated lasso based model paths for regression, logistic and multinomial regression using difference convex technique and coordinate descent. The algorithm is extremely fast, and exploits sparsity in the input x matrix where it exists. A variety of predictions can be made from the fitted models.
}
\details{
\tabular{ll}{
Package: \tab glmtlp\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2018-02-01\cr
License: \tab GPL-2 \cr
}
Very simple to use. Accepts \code{x,y} data for regression models, and
produces the regularization path over a grid of values for the tuning
parameter \code{lambda} and \code{tau}. Only 5 functions:
\code{glmtlp}\cr
\code{predict.glmnet}\cr
\code{plot.glmnet}\cr
\code{print.glmnet}\cr
\code{coef.glmnet}
}
\author{Chong Wu, Wei Pan\cr
Maintainer: Chong Wu <wuxx0845@umn.edu>
}
\references{Xiaotong Shen , Wei Pan and Yunzhang Zhu (2012)
  \emph{Likelihood-Based Selection and Sharp Parameter Estimation}, 
  \emph{Journal of the American Statistical Association, 107:497, 223-232}\cr
}
\keyword{models}
\keyword{regression}
\keyword{ package }
