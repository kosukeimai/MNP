\name{summary.mnp}

\alias{summary.mnp}
\alias{print.summary.mnp}

\title{Summarizing the results for the Multinomial Probit Models}

\description{
  \code{summary} method for class \code{mnp}.
}

\usage{
  \method{summary}{mnp}(object, CI=c(2.5, 97.5), ...)

  \method{print}{summary.mnp}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{An output object from \code{mnp}.}
  \item{CI}{A 2 dimensional vector of lower and upper bounds for the
    credible intervals used to summarize the results. The default is the
    equal tail 95 percent credible interval.
  }
  \item{x}{An object of class \code{summary.mnp}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  \code{summary.mnp} yields an object of class \code{summary.mnp}
  containing the following elements:
  \item{call}{The call from \code{mnp}.}
  \item{n.alt}{The total number of alternatives.}
  \item{base}{The base category used for fitting.}
  \item{n.obs}{The number of observations.}
  \item{n.param}{The number of estimated parameters.}
  \item{n.draws}{The number of Gibbs draws used for the summary.}
  \item{coef.table}{The summary of the posterior distribution of the
    coefficients. }
  \item{cov.table}{The summary of the posterior distribution of the
    covariance matrix.}
  This object can be printed by \code{print.summary.mnp}
}

\seealso{\code{mnp}; MNP home page at
  \url{http://imai.princeton.edu/research/MNP.html}}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}
}

\keyword{methods}
