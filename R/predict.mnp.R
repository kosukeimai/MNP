#' Posterior Prediction under the Bayesian Multinomial Probit Models
#' 
#' Obtains posterior predictions under a fitted (Bayesian) multinomial probit
#' model. \code{predict} method for class \code{mnp}.
#' 
#' The posterior predictive values are computed using the Monte Carlo sample
#' stored in the \code{mnp} output (or other sample if \code{newdraw} is
#' specified). Given each Monte Carlo sample of the parameters and each vector
#' of predictor variables, we sample the vector-valued latent variable from the
#' appropriate multivariate Normal distribution. Then, using the sampled
#' predictive values of the latent variable, we construct the most preferred
#' choice as well as the ordered preferences. Averaging over the Monte Carlo
#' sample of the preferred choice, we obtain the predictive probabilities of
#' each choice being most preferred given the values of the predictor
#' variables.  Since the predictive values are computed via Monte Carlo
#' simulations, each run may produce somewhat different values. The computation
#' may be slow if predictions with many values of the predictor variables are
#' required and/or if a large Monte Carlo sample of the model parameters is
#' used. In either case, setting \code{verbose = TRUE} may be helpful in
#' monitoring the progress of the code.
#' 
#' @param object An output object from \code{mnp}.
#' @param newdata An optional data frame containing the values of the predictor
#' variables. Predictions for multiple values of the predictor variables can be
#' made simultaneously if \code{newdata} has multiple rows. The default is the
#' original data frame used for fitting the model.
#' @param newdraw An optional matrix of MCMC draws to be used for posterior
#' predictions. The default is the original MCMC draws stored in \code{object}.
#' @param n.draws The number of additional Monte Carlo draws given each MCMC
#' draw of coefficients and covariance matrix. The specified number of latent
#' variables will be sampled from the multivariate normal distribution, and the
#' quantities of interest will be calculated by averaging over these draws.
#' This will be particularly useful calculating the uncertainty of predicted
#' probabilities. The default is \code{1}.
#' @param type The type of posterior predictions required. There are four
#' options: \code{type = "prob"} returns the predictive probabilities of being
#' the most preferred choice among the choice set.  \code{type = "choice"}
#' returns the Monte Carlo sample of the most preferred choice, and \code{type
#' = "order"} returns the Monte Carlo sample of the ordered preferences,
#' @param verbose logical. If \code{TRUE}, helpful messages along with a
#' progress report on the Monte Carlo sampling from the posterior predictive
#' distributions are printed on the screen. The default is \code{FALSE}.
#' @param ... additional arguments passed to other methods.
#' @return \code{predict.mnp} yields a list of class 
#' \code{predict.mnp} containing at least one of the following elements: 
#' \item{o}{A three dimensional array of the Monte Carlo sample from the posterior predictive
#' distribution of the ordered preferences. The first dimension corresponds to
#' the rows of \code{newdata} (or the original data set if \code{newdata} is
#' left unspecified), the second dimension corresponds to the alternatives in
#' the choice set, and the third dimension indexes the Monte Carlo sample. If
#' \code{n.draws} is greater than 1, then each entry will be an average over
#' these additional draws.  } 
#' \item{p}{A two or three dimensional array of the
#' posterior predictive probabilities for each alternative in the choice set
#' being most preferred. The first demension corresponds to the rows of
#' \code{newdata} (or the original data set if \code{newdata} is left
#' unspecified), the second dimension corresponds to the alternatives in the
#' choice set, and the third dimension (if it exists) indexes the Monte Carlo
#' sample. If \code{n.draws} is greater than 1, then the third diemsion exists
#' and indexes the Monte Carlo sample.  } 
#' \item{y}{A matrix of the Monte Carlo
#' sample from the posterior predictive distribution of the most preferred
#' choice. The first dimension correspond to the rows of \code{newdata} (or the
#' original data set if \code{newdata} is left unspecified), and the second
#' dimension indexes the Monte Carlo sample. \code{n.draws} will be set to 1
#' when computing this quantity of interest.  } 
#' \item{x}{A matrix of covariates
#' used for prediction }
#' @author Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@Princeton.Edu}
#' @seealso \code{mnp}
#' @keywords methods
#' @method predict mnp
#' @exportS3Method predict mnp
predict.mnp <- function(object, newdata = NULL, newdraw = NULL, n.draws = 1,
                        type = c("prob", "choice", "order"),
                        verbose = FALSE, ...){

  if (NA %in% match(type, c("prob", "choice", "order")))
    stop("Invalid input for `type'.")
  if (n.draws < 1)
    stop("Invalid input for `n.draws'.")

  p <- object$n.alt
  if (is.null(newdraw)) 
    param <- object$param
  else
    param <- newdraw
  n.cov <- ncol(coef(object))
  coef <- param[,1:n.cov]
  n.mcmc <- nrow(coef)
  cov <- param[,(n.cov+1):ncol(param)]
  
  ## get X matrix ready
  if (is.null(newdata)) 
    x <- object$x
  else {
    call <- object$call
    x <- xmatrix.mnp(stats::as.formula(call$formula), data = newdata,
                     choiceX = call$choiceX,
                     cXnames = eval(call$cXnames),
                     base = object$base, n.dim = p-1,
                     lev = object$alt, MoP = is.matrix(object$y),
                     verbose = FALSE, extra = FALSE)
    if (nrow(x) > 1) 
      x <- as.matrix(x[apply(is.na(x), 1, sum)==0,])
    else if (sum(is.na(x))>0)
      stop("Invalid input for `newdata'.")
  }
  n.obs <- nrow(x)
  if (verbose) {
    if (n.obs == 1)
      cat("There is one observation to predict. Please wait...\n")
    else
      cat("There are", n.obs, "observations to predict. Please wait...\n")
  }

  alt <- object$alt
  if (object$base != alt[1]) 
    alt <- c(object$base, alt[1:(length(alt)-1)])

  res <- .C("predict", as.double(x), as.integer(n.obs),
            as.double(coef), as.double(cov), as.integer(p-1),
            as.integer(n.cov), as.integer(n.mcmc),
            as.integer(n.draws), as.integer(verbose),
            prob = if (n.draws > 1) double(n.obs*n.mcmc*p)
                   else double(n.obs*p),
            choice = double(n.obs*n.mcmc),
            order = double(n.obs*n.mcmc*p),
            PACKAGE = "MNP")

  ans <- list()
  if ("choice" %in% type)
    ans$y <- matrix(res$choice, ncol=n.mcmc, nrow = n.obs,
                    byrow = TRUE, dimnames = list(rownames(x), 1:n.mcmc))
  else
    ans$y <- NULL
  if ("order" %in% type)
    ans$o <- aperm(array(res$order, dim=c(p, n.mcmc, n.obs),
                         dimnames = list(alt, 1:n.mcmc,
                           rownames(x))), c(3,1,2))
  else
    ans$o <- NULL
  if ("prob" %in% type)
    if (n.draws > 1)
      ans$p <- aperm(array(res$prob, dim = c(p, n.mcmc, n.obs),
                           dimnames = list(alt, 1:n.mcmc,
                             rownames(x))), c(3,1,2))
    else
      ans$p <- matrix(res$prob, nrow = n.obs, ncol = p, byrow = TRUE,
                      dimnames = list(rownames(x), alt))
  else
    ans$p <- NULL

  ans$x <- x
  class(ans) <- "predict.mnp"
  return(ans)
}
