#' Extract Multinomial Probit Model Covariance Matrix
#'
#' \code{vcov.mnp} is a function which extracts the posterior draws of
#' covariance matrix from objects returned by \code{mnp}.
#'
#'
#' @param object An output object from \code{mnp}.
#' @param subset A scalar or a numerical vector specifying the row number(s) of
#' \code{param} in the output object from \code{mnp}. If specified, the
#' posterior draws of covariance matrix for those rows are extracted. The
#' default is \code{NULL} where all the posterior draws are extracted.
#' @param ... further arguments passed to or from other methods.
#' @return When a numerical vector or \code{NULL} is specified for
#' \code{subset} argument, \code{vcov.mnp} returns a three dimensional array
#' where the third dimension indexes posterior draws. When a scalar is
#' specified for \code{subset} arugment, \code{vcov.mnp} returns a matrix.
#' @author Kosuke Imai, Department of Government, Harvard University
#' \email{imai@@harvard.edu}
#' @seealso \code{mnp}, \code{coef.mnp};
#' @keywords methods
#' @exportS3Method vcov mnp
vcov.mnp <- function(object, subset = NULL, ...) {
  if (is.null(subset))
    subset <- 1:nrow(object$param)
  else if (max(subset) > nrow(object$param))
    stop(paste("invalid input for `subset.' only", nrow(object$param), "draws are stored."))

  p <- object$n.alt
  n <- length(subset)
  Sigma <- array(0, c(p-1, p-1, n))
  param <- object$param
  n.cov <- ncol(param) - p*(p-1)/2
  cov <- param[,(n.cov+1):ncol(param)]
  for (i in 1:n) {
    count <- 1
    for (j in 1:(p-1)) {
      Sigma[j,j:(p-1),i] <- cov[subset[i],count:(count+p-j-1)]
      count <- count + p - j
    }
    diag(Sigma[,,i]) <- diag(Sigma[,,i]/2)
    Sigma[,,i] <- Sigma[,,i] + t(Sigma[,,i])
  }
  tmp <- list()
  tmp[[1]] <- tmp[[2]] <- object$alt[-pmatch(object$base, object$alt)]
  tmp[[3]] <- as.character(subset)
  dimnames(Sigma) <- tmp
  if (n > 1)
    return(Sigma)
  else
    return(Sigma[,,1])
}
