\name{cvamEstimate}
\alias{cvamEstimate}
\alias{print.cvamEstimate}
\alias{print.cvamEstimateList}

\title{Obtain Estimated Probabilities from  a Fitted Model}

\description{
After fitting a log-linear model with \code{cvam}, 
the fitted model object may be passed to this function to obtain estimated
marginal and conditional probabilities for model factors.
}

\usage{
cvamEstimate(estimate, obj, meanSeries = TRUE,
   confidence = obj$control$confidence,
   probRound = obj$control$probRound, probDigits = obj$control$probDigits, ...)

\method{print}{cvamEstimate}(x, showHeader = TRUE, ...)

\method{print}{cvamEstimateList}(x, showHeader = TRUE, ...)
}

\arguments{
  \item{estimate}{a formula or list of formulas indicating the desired
     probabilities; see DETAILS.}

  \item{obj}{an object produced by \code{\link{cvam}} containing 
  results from a model fit}

  \item{meanSeries}{applies when \code{obj} contains results from a
     simulation run. 
     If \code{TRUE}, then the requested estimates are computed based on
     a running mean of cell probabilities over all iterations 
     after the burn-in period. If \code{FALSE}, then the requested
     estimates are based only on the cell probabilities from the final
     iteration, and (assuming the run was sufficiently long, if it is MCMC)
     can be regarded as a single draw from their posterior distribution.} 

  \item{confidence}{confidence coefficient for asymmetric interval
     estimates; see DETAILS.}

  \item{probRound}{if TRUE, probabilities will be rounded.}

  \item{probDigits}{number of digits for rounding probabilities.}

  \item{x}{a set of estimates to be printed.}

  \item{showHeader}{if \code{TRUE}, a descriptive header is printed.}

  \item{...}{additional arguments to be passed to \code{print}.}
}

\details{
The argument \code{estimate} should be a one-sided formula or a list
of one-sided formulas, with variables separated by `\code{+}', and
variables to be conditioned on appearing after `\code{|}'. For
example, \code{~ A} requests marginal probabilities for every level of
\code{A};
\code{~ A + B | C + D} requests conditional probabilities for
every level combination of \code{A} and \code{B} given every level
combination of \code{C} and \code{D}.

\itemize{
\item If \code{obj} was produced with \code{saturated=FALSE} and
\code{method="EM"}, then standard errors for all 
probabilities are computed using Taylor linearization, also known as
the delta method, based on the asymptotic covariance matrix for the
log-linear coefficients.
\item If \code{obj} was produced with \code{saturated=FALSE} and
\code{method="MCMC"} or \code{"approxBayes"}, then standard errors are
computed with Taylor 
linearization, based on the covariance matrix for the simulated log-linear
coefficients from all iterations after the burn-in period.
\item If \code{obj} was produced with \code{saturated=TRUE}, then 
standard errors are not computed.
}

A symmetric confidence interval for a
probability may be problematic, especially if the estimate is close to
zero or one. Asymmetric confidence intervals are computed by applying
a normal approximation on the logistic (log-odds) scale and
translating the endpoints back to the probability scale.
}

\value{
if \code{estimate} is a single formula, this function returns a data
frame containing 
estimated probabilities, standard errors, and endpoints of approximate
confidence intervals. If \code{estimate} is a list of formulas, then a
list of data frames is returned.
}

\note{
Estimates may also be requested when fitting a model with \code{cvam},
by providing a formula or list of formulas to the optional argument
\code{estimate}.
}

\references{
For more information, refer to the package vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}}
\code{\link{cvamPredict}}
\code{\link{cvamImpute}}
\code{\link{cvamLik}}
}

\examples{
fit <- cvam( ~ Sex * PolViews * AbAny, data=abortion2000 )
cvamEstimate( list( ~ AbAny | Sex, ~ AbAny | PolViews ), fit )
}
