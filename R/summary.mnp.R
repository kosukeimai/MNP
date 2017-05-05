#' Summarizing the results for the Multinomial Probit Models
#' 
#' \code{summary} method for class \code{mnp}.
#' 
#' 
#' @aliases summary.mnp print.summary.mnp
#' @param object An output object from \code{mnp}.
#' @param CI A 2 dimensional vector of lower and upper bounds for the credible
#' intervals used to summarize the results. The default is the equal tail 95
#' percent credible interval.
#' @param x An object of class \code{summary.mnp}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#' @return \code{summary.mnp} yields an object of class \code{summary.mnp}
#' containing the following elements: \item{call}{The call from \code{mnp}.}
#' \item{n.alt}{The total number of alternatives.} \item{base}{The base
#' category used for fitting.} \item{n.obs}{The number of observations.}
#' \item{n.param}{The number of estimated parameters.} \item{n.draws}{The
#' number of Gibbs draws used for the summary.} \item{coef.table}{The summary
#' of the posterior distribution of the coefficients. } \item{cov.table}{The
#' summary of the posterior distribution of the covariance matrix.} This object
#' can be printed by \code{print.summary.mnp}
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@Princeton.Edu}
#' @seealso \code{mnp}; MNP home page at
#' \url{http://imai.princeton.edu/research/MNP.html}
#' @keywords methods
#' @export summary.mnp
summary.mnp <- function(object, CI=c(2.5, 97.5),...){

  p <- object$n.alt
  param <- object$param
  n.cov <- ncol(param) - p*(p-1)/2
  n.draws <- nrow(param)
  param.table <- cbind(apply(param, 2, mean), apply(param, 2, sd),
                       apply(param, 2, quantile, min(CI)/100),
                       apply(param, 2, quantile, max(CI)/100)) 
  colnames(param.table) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                             paste(max(CI), "%", sep=""))
  rownames(param.table) <- colnames(param)
  
  ans <- list(call = object$call, base = object$base, n.alt=p,
              n.obs = if(is.matrix(object$y)) nrow(object$y) else length(object$y),
              n.param = ncol(param)-1, n.draws = n.draws,
              coef.table= if(n.cov > 1) param.table[1:n.cov,]
              else matrix(param.table[1,], nrow=1,
                          dimnames=list(rownames(param.table)[1], colnames(param.table))),
              cov.table=param.table[(n.cov+1):ncol(param),])  
  class(ans) <- "summary.mnp"
  return(ans)
}
