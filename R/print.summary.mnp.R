#' Print the summary of the results for the Multinomial Probit Models
#'
#' \code{summary} print method for class \code{mnp}.
#'
#' @param x An object of class \code{summary.mnp}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@Princeton.Edu}
#' @seealso \code{mnp}
#' @keywords methods
#' @exportS3Method print summary.mnp
print.summary.mnp <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      sep = "")

  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coef.table, digits = digits, na.print = "NA", ...)

  cat("\nCovariances:\n")
  stats::printCoefmat(x$cov.table, digits = digits, na.print = "NA", ...)

  cat("\nBase category:", x$base)
  cat("\nNumber of alternatives:", x$n.alt)
  cat("\nNumber of observations:", x$n.obs)
  cat("\nNumber of estimated parameters:", x$n.param)
  cat("\nNumber of stored MCMC draws:", x$n.draws)
  cat("\n\n")
  invisible(x)
}
