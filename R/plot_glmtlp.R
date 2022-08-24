#' Plot Method for a "glmtlp" Object
#'
#' @description
#' Generates a solution path plot for a fitted \code{"glmtlp"} object.
#'
#' @details
#' The generated plot is a \code{ggplot} object, and therefore, the users are able
#'   to customize the plots following the \code{ggplot2} syntax.
#'
#' @param x Fitted \code{glmtlp} object.
#' @param xvar The x-axis variable to plot against, including \code{"lambda"},
#'   \code{"kappa"}, \code{"deviance"}, \code{"l1_norm"}, and \code{"log_lambda"}.
#' @param xlab The x-axis label of the plot, default is \code{"Lambda"},
#'   \code{"Kappa"}, \code{"Fraction of Explained Deviance"}, \code{"L1 Norm"},
#'   and \code{"Log Lambda"}.
#' @param ylab The y-axis label of the plot, default is "Coefficients".
#' @param title The main title of the plot, default is "Solution Path".
#' @param label Logical, whether or not attach the labels for the non-zero
#'   coefficients, default is \code{FALSE}.
#' @param label.size The text size of the labels, default is 3.
#' @param \dots Additional arguments.
#'
#' @return A \code{ggplot} object.
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
#'   \cr Shen, X., Pan, W., Zhu, Y., & Zhou, H. (2013).
#'   \emph{On constrained and regularized high-dimensional regression.
#'   Annals of the Institute of Statistical Mathematics, 65(5), 807-832.}
#'   \cr Li, C., Shen, X., & Pan, W. (2021).
#'   \emph{Inference for a Large Directed Graphical Model with Interventions.
#'   arXiv preprint arXiv:2110.03805.}
#'   \cr Yang, Y., & Zou, H. (2014).
#'   \emph{A coordinate majorization descent algorithm for l1 penalized learning.
#'   Journal of Statistical Computation and Simulation, 84(1), 84-95.}
#'   \cr Two R package Github: \emph{ncvreg} and \emph{glmnet}.
#'
#' @keywords models plot
#'
#' @examples
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' fit <- glmtlp(X, y, family = "gaussian", penalty = "l1")
#' plot(fit, xvar = "lambda")
#' plot(fit, xvar = "log_lambda")
#' plot(fit, xvar = "l1_norm")
#' plot(fit, xvar = "log_lambda", label = TRUE)
#' fit2 <- glmtlp(X, y, family = "gaussian", penalty = "l0")
#' plot(fit2, xvar = "kappa", label = TRUE)
#'
#' @import ggplot2
#' @method plot glmtlp
#' @export
#' @export plot.glmtlp

plot.glmtlp <- function(x,xvar=c("lambda", "kappa", "deviance", "l1_norm", "log_lambda"),
                        xlab=iname, ylab="Coefficients", title="Solution Path",
                        label=FALSE, label_size=3, ...) {
    xvar <- match.arg(xvar)
    which_plot <- which(apply(abs(x$beta), 1, sum) != 0)
    if (length(which_plot) == 0) {
        warning("No plot produced, since all coefficients are zero or non-penalized.")
        return ()
    }
    if (xvar == "lambda" || xvar == "log_lambda") {
        if (x$method == "tlp-constrained") stop("No plot generated, since the l0 penalty should be plotted with xvar='kappa'.")
        if (length(x$lambda) == 1) {
        warning("No plot generated, since the x was fit only on one lambda.")
        return ()
        }
    }
    if (xvar == "kappa") {
        if (x$method != "tlp-constrained") stop("No plot generated, since the xvar='kappa' should be used together with non-l0 penalties.")
        if (length(x$kappa) == 1) {
        warning("No plot generated, since the x was fit only on one kappa.")
        return ()
        }
    }

    beta <- as.matrix(x$beta)[which_plot, , drop=FALSE]
    #xvar <- match.arg(xvar)
    switch(xvar,
        "lambda" = {
            index <- x$lambda
            iname <- expression(lambda)
        },
        "kappa" = {
            index <- x$kappa
            iname <- expression(kappa)
        },
        "deviance" = {
            index <- x$deviance
            iname <- "Fraction of Explained Deviance"
        },
        "l1_norm" = {
            index <- apply(abs(beta), 2, sum)
            iname <- "L1 Norm"
        },
        "log_lambda" = {
            index <- log(x$lambda)
            iname <- expression(Log(lambda))
        }
    )

    df <- data.frame(index = rep(index, each = nrow(beta)),
                     variable = rep(rownames(beta), ncol(beta)),
                     value = c(beta))
    df_legend <- data.frame(matrix(nrow=nrow(beta), ncol=3))
    colnames(df_legend) <- c("x_pos", "y_pos", "label")
    if (xvar == 'kappa') {
      df_legend$x_pos <- max(index)
      hjust_val <- 0
    } else {
      df_legend$x_pos <- min(index)
      hjust_val <- 1
    }
    df_legend$y_pos <- beta[, ncol(beta)]
    df_legend$label <- rownames(beta)
    g <- ggplot(df, aes_string(x = "index", y = "value")) +
      geom_line(aes_string(colour = "variable", group = "variable")) +
      theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
      xlab(xlab) + ylab(ylab) + ggtitle(title)
    if (label) {
      g <- g + geom_text(data=df_legend, aes_string(x = "x_pos", y = "y_pos", label = "label", colour = "label",
                                             hjust = hjust_val), size = label_size)
    }
    if (xvar == "kappa") {
      g <- g + scale_x_discrete(limits = factor(x$kappa))
    }

    g
}
