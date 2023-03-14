#' Extract Multinomial Probit Model Coefficients
#'
#' \code{coef.mnp} is a function which extracts multinomial probit model
#' coefficients from ojbects returned by \code{mnp}. \code{coefficients.mnp} is
#' an alias for it. \code{coef} method for class \code{mnp}.
#'
#'
#' @param object An output object from \code{mnp}.
#' @param subset A scalar or a numerical vector specifying the row number(s) of
#' \code{param} in the output object from \code{mnp}. If specified, the
#' posterior draws of coefficients for those rows are extracted. The default is
#' \code{NULL} where all the posterior draws are extracted.
#' @param ... further arguments passed to or from other methods.
#' @return \code{coef.mnp} returns a matrix (when a numerical vector or
#' \code{NULL} is specified for \code{subset} argument) or a vector (when a
#' scalar is specified for \code{subset} arugment) of multinomila probit model
#' coefficients.
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@Princeton.Edu}
#' @seealso \code{mnp}, \code{vcov.mnp};
#' @keywords methods
#' @exportS3Method coef mnp
coef.mnp <- function(object, subset = NULL, ...) {
  param <- object$param
  p <- object$n.alt
  n.cov <- ncol(param) - p*(p-1)/2
  n <- nrow(param)
  if (is.null(subset))
    return(param[,1:n.cov])
  else if (subset > n)
    stop(paste("invalid input for `subset.' only", nrow(param), "draws are stored."))
  else
    return(param[subset, 1:n.cov])
}
