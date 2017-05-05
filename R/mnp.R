#' Fitting the Multinomial Probit Model via Markov chain Monte Carlo
#' 
#' \code{mnp} is used to fit (Bayesian) multinomial probit model via Markov
#' chain Monte Carlo.  \code{mnp} can also fit the model with different choice
#' sets for each observation, and complete or partial ordering of all the
#' available alternatives. The computation uses the efficient marginal data
#' augmentation algorithm that is developed by Imai and van Dyk (2005a).
#' 
#' To fit the multinomial probit model when only the most preferred choice is
#' observed, use the syntax for the formula, \code{y ~ x1 + x2}, where \code{y}
#' is a factor variable indicating the most preferred choice and \code{x1} and
#' \code{x2} are individual-specific covariates. The interactions of
#' individual-specific variables with each of the choice indicator variables
#' will be fit.
#' 
#' To specify choice-specific covariates, use the syntax,
#' \code{choiceX=list(A=cbind(z1, z2), B=cbind(z3, z4), C=cbind(z5, z6))},
#' where \code{A}, \code{B}, and \code{C} represent the choice names of the
#' response variable, and \code{z1} and \code{z2} are each vectors of length
#' \eqn{n} that record the values of the two choice-specific covariates for
#' each individual for choice A, likewise for \code{z3}, \eqn{\ldots},
#' \code{z6}. The corresponding variable names via \code{cXnames=c("price",
#' "quantity")} need to be specified, where \code{price} refers to the
#' coefficient name for \code{z1}, \code{z3}, and \code{z5}, and
#' \code{quantity} refers to that for \code{z2}, \code{z4}, and \code{z6}.
#' 
#' If the choice set varies from one observation to another, use the syntax,
#' \code{cbind(y1, y2, y3) ~ x1 + x2}, in the case of a three choice problem,
#' and indicate unavailable alternatives by \code{NA}. If only the most
#' preferred choice is observed, \code{y1}, \code{y2}, and \code{y3} are
#' indicator variables that take on the value one for individuals who prefer
#' that choice and zero otherwise. The last column of the response matrix,
#' \code{y3} in this particular example syntax, is used as the base category.
#' 
#' To fit the multinomial probit model when the complete or partial ordering of
#' the available alternatives is recorded, use the same syntax as when the
#' choice set varies (i.e., \code{cbind(y1, y2, y3, y4) ~ x1 + x2}). For each
#' observation, all the available alternatives in the response variables should
#' be numerically ordered in terms of preferences such as \code{1 2 2 3}. Ties
#' are allowed. The missing values in the response variable should be denoted
#' by \code{NA}. The software will impute these missing values using the
#' specified covariates. The resulting uncertainty estimates of the parameters
#' will properly reflect the amount of missing data. For example, we expect the
#' standard errors to be larger when there is more missing data.
#' 
#' @aliases mnp MNP
#' @param formula A symbolic description of the model to be fit specifying the
#' response variable and covariates. The formula should not include the
#' choice-specific covariates. Details and specific examples are given below.
#' @param data An optional data frame in which to interpret the variables in
#' \code{formula} and \code{choiceX}. The default is the environment in which
#' \code{mnp} is called.
#' @param choiceX An optional list containing a matrix of choice-specific
#' covariates for each category. Details and examples are provided below.
#' @param cXnames A vector of the names for the choice-specific covariates
#' specified in \code{choiceX}. The details and examples are provided below.
#' @param base The name of the base category. For the standard multinomial
#' probit model, the default is the lowest level of the response variable. For
#' the multinomial probit model with ordered preferences, the default base
#' category is the last column in the matrix of response variables.
#' @param latent logical. If \code{TRUE}, then the latent variable W will be
#' returned. See Imai and van Dyk (2005) for the notation. The default is
#' \code{FALSE}.
#' @param invcdf logical. If \code{TRUE}, then the inverse cdf method is used
#' for truncated normal sampling. If \code{FALSE}, then the rejection sampling
#' method is used. The default is \code{FALSE}.
#' @param trace logical. If \code{TRUE}, then the trace of the variance
#' covariance matrix is set to a constant (here, it is equal to \code{n.dim})
#' instead of setting its first diagonal element to 1.  The former avoids the
#' arbitrariness of fixing one particular diagonal element in order to achieve
#' identification (see Burgette and Nordheim, 2009).
#' @param n.draws A positive integer. The number of MCMC draws. The default is
#' \code{5000}.
#' @param p.var A positive definite matrix. The prior variance of the
#' coefficients.  A scalar input can set the prior variance to the diagonal
#' matrix whose diagonal element is equal to that value. The default is
#' \code{"Inf"}, which represents an improper noninformative prior distribution
#' on the coefficients.
#' @param p.df A positive integer greater than \code{n.dim-1}. The prior
#' degrees of freedom parameter for the covariance matrix. The default is
#' \code{n.dim+1}, which is equal to the total number of alternatives.
#' @param p.scale A positive definite matrix.  When \code{trace = FALSE}, its
#' first diagonal element is set to \code{1} if it is not equal to 1 already.
#' The prior scale matrix for the covariance matrix. A scalar input can be used
#' to set the scale matrix to a diagonal matrix with diagonal elements equal to
#' the scalar input value. The default is \code{1}.
#' @param coef.start A vector. The starting values for the coefficients.  A
#' scalar input sets the starting values for all the coefficients equal to that
#' value.  The default is \code{0}.
#' @param cov.start A positive definite matrix. When \code{trace = FALSE}, its
#' first diagonal element is set to \code{1} if it is not equal to 1 already.
#' The starting values for the covariance matrix. A scalar input can be used to
#' set the starting value to a diagonal matrix with diagonal elements equal to
#' the scalar input value. The default is \code{1}.
#' @param burnin A positive integer. The burnin interval for the Markov chain;
#' i.e., the number of initial Gibbs draws that should not be stored. The
#' default is \code{0}.
#' @param thin A positive integer. The thinning interval for the Markov chain;
#' i.e., the number of Gibbs draws between the recorded values that are
#' skipped. The default is \code{0}.
#' @param verbose logical. If \code{TRUE}, helpful messages along with a
#' progress report of the Gibbs sampling are printed on the screen. The default
#' is \code{FALSE}.
#' @return An object of class \code{mnp} containing the following elements:
#' \item{param}{A matrix of the Gibbs draws for each parameter; i.e., the
#' coefficients and covariance matrix. For the covariance matrix, the elements
#' on or above the diagonal are returned.  } 
#' \item{call}{The matched call.}
#' \item{x}{The matrix of covariates.} 
#' \item{y}{The vector or matrix of the
#' response variable.} 
#' \item{w}{The three dimensional array of the latent
#' variable, W. The first dimension represents the alternatives, and the second
#' dimension indexes the observations. The third dimension represents the Gibbs
#' draws. Note that the latent variable for the base category is set to 0, and
#' therefore omitted from the output.} 
#' \item{alt}{The names of alternatives.}
#' \item{n.alt}{The total number of alternatives.} 
#' \item{base}{The base
#' category used for fitting.} 
#' \item{invcdf}{The value of 
#' \code{invcdf}.}
#' \item{p.var}{The prior variance for the coefficients.} 
#' \item{p.df}{The prior
#' degrees of freedom parameter for the covariance matrix.} 
#' \item{p.scale}{The
#' prior scale matrix for the covariance matrix.} 
#' \item{burnin}{The number of
#' initial burnin draws.} 
#' \item{thin}{The thinning interval.}
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@Princeton.Edu}, \url{http://imai.princeton.edu}; David A. van
#' Dyk, Department of Statistics, University of California, Irvine
#' \email{dvd@@uci.edu}, \url{http://www.ics.uci.edu/~dvd}.
#' @seealso \code{coef.mnp}, \code{cov.mnp}, \code{predict.mnp},
#' \code{summary.mnp}; MNP home page at
#' \url{http://imai.princeton.edu/research/MNP.html}
#' @references Imai, Kosuke and David A. van Dyk. (2005a) \dQuote{A Bayesian
#' Analysis of the Multinomial Probit Model Using the Marginal Data
#' Augmentation,} \emph{Journal of Econometrics}, Vol. 124, No. 2 (February),
#' pp.311-334.
#' 
#' Imai, Kosuke and David A. van Dyk. (2005b) \dQuote{MNP: R Package for
#' Fitting the Multinomial Probit Models,} \emph{Journal of Statistical
#' Software}, Vol. 14, No. 3 (May), pp.1-32.
#' 
#' Burgette, L.F. and E.V. Nordheim. (2009).  \dQuote{An alternate identifying
#' restriction for the Bayesian multinomial probit model,} \emph{Technical
#' report}, Department of Statistics, University of Wisconsin, Madison.
#' @keywords models
#' @examples
#' 
#' ###
#' ### NOTE: this example is not fully analyzed. In particular, the
#' ### convergence has not been assessed. A full analysis of these data
#' ### sets appear in Imai and van Dyk (2005b).
#' ###
#' 
#' ## load the detergent data
#' data(detergent)
#' ## run the standard multinomial probit model with intercepts and the price
#' res1 <- mnp(choice ~ 1, choiceX = list(Surf=SurfPrice, Tide=TidePrice,
#'                                        Wisk=WiskPrice, EraPlus=EraPlusPrice,
#'                                        Solo=SoloPrice, All=AllPrice),
#'             cXnames = "price", data = detergent, n.draws = 500, burnin = 100,
#'             thin = 3, verbose = TRUE)
#' ## summarize the results
#' summary(res1)
#' ## calculate the quantities of interest for the first 3 observations
#' pre1 <- predict(res1, newdata = detergent[1:3,])
#' 
#' ## load the Japanese election data
#' data(japan)
#' ## run the multinomial probit model with ordered preferences
#' res2 <- mnp(cbind(LDP, NFP, SKG, JCP) ~ gender + education + age, data = japan,
#'             verbose = TRUE)
#' ## summarize the results
#' summary(res2)
#' ## calculate the predicted probabilities for the 10th observation
#' ## averaging over 100 additional Monte Carlo draws given each of MCMC draw.
#' pre2 <- predict(res2, newdata = japan[10,], type = "prob", n.draws = 100,
#'                 verbose = TRUE)
#' 
#' @export mnp
mnp <- function(formula, data = parent.frame(), choiceX = NULL,
                cXnames = NULL, base = NULL, latent = FALSE,
                invcdf = FALSE, trace = TRUE, n.draws = 5000, p.var = "Inf", 
                p.df = n.dim+1, p.scale = 1, coef.start = 0,
                cov.start = 1, burnin = 0, thin = 0, verbose = FALSE) {   
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$choiceX <- mf$cXnames <- mf$base <- mf$n.draws <- mf$latent <-
    mf$p.var <- mf$p.df <- mf$p.scale <- mf$coef.start <- mf$invcdf <-
      mf$trace <- mf$cov.start <- mf$verbose <- mf$burnin <- mf$thin <- NULL   
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)

  ## fix this parameter
  p.alpha0 <- 1
  
  ## obtaining Y
  tmp <- ymatrix.mnp(mf, base=base, extra=TRUE, verbose=verbose)
  Y <- tmp$Y
  MoP <- tmp$MoP
  lev <- tmp$lev
  base <- tmp$base
  p <- tmp$p
  n.dim <- p - 1
  if(verbose)
    cat("\nThe base category is `", base, "'.\n\n", sep="") 
  if (p < 3)
    stop("The number of alternatives should be at least 3.")
  if(verbose) 
    cat("The total number of alternatives is ", p, ".\n\n", sep="") 
  if(verbose) {
    if (trace)
      cat("The trace restriction is used instead of the diagonal restriction.\n\n")
    else
      cat("The diagonal restriction is used instead of the trace restriction.\n\n")
  }
  
  ### obtaining X
  tmp <- xmatrix.mnp(formula, data=eval.parent(data),
                     choiceX=call$choiceX, cXnames=cXnames, 
                     base=base, n.dim=n.dim, lev=lev, MoP=MoP,
                     verbose=verbose, extra=TRUE)
  X <- tmp$X
  coefnames <- tmp$coefnames
  n.cov <- ncol(X) / n.dim
  
  ## listwise deletion for X
  na.ind <- apply(is.na(X), 1, sum)
  if (ncol(Y) == 1)
    na.ind <- na.ind + is.na(Y)
  Y <- Y[na.ind==0,]
  X <- X[na.ind==0,]
  n.obs <- nrow(X)
  
  if (verbose) {
    cat("The dimension of beta is ", n.cov, ".\n\n", sep="")
    cat("The number of observations is ", n.obs, ".\n\n", sep="")
    if (sum(na.ind>0)>0) {
      if (sum(na.ind>0)==1)
        cat("The observation ", (1:length(na.ind))[na.ind>0], " is dropped due to missing values.\n\n", sep="")
      else {
        cat("The following ", sum(na.ind>0), " observations are dropped due to missing values:\n", sep="")
        cat((1:length(na.ind))[na.ind>0], "\n\n")
      }
    }
  } 
  
  ## checking the prior for beta
  p.imp <- FALSE 
  if (p.var == Inf) {
    p.imp <- TRUE
    p.prec <- diag(0, n.cov)
    if (verbose)
      cat("Improper prior will be used for beta.\n\n")
  }
  else if (is.matrix(p.var)) {
    if (ncol(p.var) != n.cov || nrow(p.var) != n.cov)
      stop("The dimension of `p.var' should be ", n.cov, " x ", n.cov, sep="")
    if (sum(sign(eigen(p.var)$values) < 1) > 0)
      stop("`p.var' must be positive definite.")
    p.prec <- solve(p.var)
  }
  else {
    p.var <- diag(p.var, n.cov)
    p.prec <- solve(p.var)
  }
  p.mean <- rep(0, n.cov)

  ## checking prior for Sigma
  p.df <- eval(p.df)
  if (length(p.df) > 1)
    stop("`p.df' must be a positive integer.")
  if (p.df < n.dim)
    stop(paste("`p.df' must be at least ", n.dim, ".", sep=""))
  if (abs(as.integer(p.df) - p.df) > 0)
    stop("`p.df' must be a positive integer.")
  if (!is.matrix(p.scale))  
    p.scale <- diag(p.scale, n.dim)
  if (ncol(p.scale) != n.dim || nrow(p.scale) != n.dim)
    stop("`p.scale' must be ", n.dim, " x ", n.dim, sep="")
  if (sum(sign(eigen(p.scale)$values) < 1) > 0)
    stop("`p.scale' must be positive definite.")
  else if ((trace == FALSE) & (p.scale[1,1] != 1)) {
    p.scale[1,1] <- 1
    warning("p.scale[1,1] will be set to 1.")
  }
  Signames <- NULL
  for(j in 1:n.dim)
    for(k in 1:n.dim)
      if (j<=k)
        Signames <- c(Signames, paste(if(MoP) lev[j] else lev[j+1],
                                      ":", if(MoP) lev[k] else lev[k+1], sep="")) 

  ## checking starting values
  if (length(coef.start) == 1)
    coef.start <- rep(coef.start, n.cov)
  else if (length(coef.start) != n.cov)
    stop(paste("The dimenstion of `coef.start' must be  ",
               n.cov, ".", sep=""))
  if (!is.matrix(cov.start)) {
    cov.start <- diag(n.dim)*cov.start
    if (!trace)
      cov.start[1,1] <- 1
  }
  else if (ncol(cov.start) != n.dim || nrow(cov.start) != n.dim)
    stop("The dimension of `cov.start' must be ", n.dim, " x ", n.dim, sep="")
  else if (sum(sign(eigen(cov.start)$values) < 1) > 0)
    stop("`cov.start' must be a positive definite matrix.")
  else if ((trace == FALSE) & (cov.start[1,1] != 1)) {
    cov.start[1,1] <- 1
    warning("cov.start[1,1] will be set to 1.")
  }
  
  ## checking thinnig and burnin intervals
  if (burnin < 0)
    stop("`burnin' should be a non-negative integer.") 
  if (thin < 0)
    stop("`thin' should be a non-negative integer.")
  keep <- thin + 1
  
  ## running the algorithm
  if (latent)
    n.par <- n.cov + n.dim*(n.dim+1)/2 + n.dim*n.obs
  else
    n.par <- n.cov + n.dim*(n.dim+1)/2
  if(verbose)
    cat("Starting Gibbs sampler...\n")
  # recoding NA into -1
  Y[is.na(Y)] <- -1

  param <- .C("cMNPgibbs", as.integer(n.dim),
              as.integer(n.cov), as.integer(n.obs), as.integer(n.draws),
              as.double(p.mean), as.double(p.prec), as.integer(p.df),
              as.double(p.scale*p.alpha0), as.double(X), as.integer(Y), 
              as.double(coef.start), as.double(cov.start), 
              as.integer(p.imp), as.integer(invcdf),
              as.integer(burnin), as.integer(keep), as.integer(trace),
              as.integer(verbose), as.integer(MoP), as.integer(latent),
              pdStore = double(n.par*floor((n.draws-burnin)/keep)),
              PACKAGE="MNP")$pdStore 
  param <- matrix(param, ncol = n.par,
                  nrow = floor((n.draws-burnin)/keep), byrow=TRUE)
  if (latent) {
    W <- array(as.vector(t(param[,(n.par-n.dim*n.obs+1):n.par])),
               dim = c(n.dim, n.obs, floor((n.draws-burnin)/keep)),
               dimnames = list(lev[!(lev %in% base)], rownames(Y), NULL))
    param <- param[,1:(n.par-n.dim*n.obs)]
    }
  else
    W <- NULL
  colnames(param) <- c(coefnames, Signames)
    
  ##recoding -1 back into NA
  Y[Y==-1] <- NA

  ## returning the object
  res <- list(param = param, x = X, y = Y, w = W, call = call, alt = lev,
              n.alt = p, base = base, invcdf = invcdf, trace = trace,
              p.mean = if(p.imp) NULL else p.mean, p.var = p.var, 
              p.df = p.df, p.scale = p.scale, burnin = burnin, thin = thin) 
  class(res) <- "mnp"
  return(res)
}
  


