\name{cvamControl}
\alias{cvamControl}

\title{
Control Parameters for \code{cvam}
}

\description{
The \code{\link{cvam}} function fits
log-linear models to
coarsened categorical variables. Its model-fitting
procedures are governed by parameters in a \code{cvamControl} object
created by the auxiliary function documented here. This function is
intended for internal use; the only reason to invoke this
function directly is to display the control parameters and their
default values.
}

\usage{
cvamControl( iterMaxEM = 500L, iterMaxNR = 50L,
   iterApproxBayes = 1L, imputeApproxBayes = FALSE,
   iterMCMC = 5000L, burnMCMC = 0L, thinMCMC = 1L, imputeEvery = 0L,
   saveProbSeries = FALSE,
   typeMCMC = c("DA","RWM"), tuneDA = c(10,.8,.8), tuneRWM = c(1000,.1),
   stuckLimit = 25L,
   startValDefault = c("center", "uniform"), startValJitter = 0,
   critEM = 1e-06, critNR = 1e-06, critBoundary = 1e-08, ncolMaxMM = 1000L,
   excludeAllNA = TRUE, critModelCheck=1e-08, confidence=.95, 
   probRound = TRUE, probDigits = 4L ) 
}

\arguments{
  \item{iterMaxEM}{maximum number of iterations performed when
  \code{method = "EM"}; see DETAILS.}

  \item{iterMaxNR}{maximum number of iterations of Newton-Raphson
  performed during an M-step of EM; see DETAILS.}

  \item{iterApproxBayes}{number of simulated log-linear coefficient
  vectors to be drawn from their approximate posterior distribution
  when \code{method="approxBayes"}.}

  \item{imputeApproxBayes}{if \code{TRUE} then, for each draw of the
  log-linear coefficients from their approximate posterior distribution,
  the true frequencies will be imputed.}

  \item{iterMCMC}{number of iterations of Markov chain Monte Carlo
  after the burn-in period when \code{method="MCMC"}.}

  \item{burnMCMC}{number of iterations of Markov chain Monte Carlo
  performed as a burn-in period, for which the reults are
  discarded. The total number of iterations performed is
  \code{burnMCMC+iterMCMC}.} 

  \item{thinMCMC}{thinning interval for saving the results from MCMC
  as a series.} 

  \item{imputeEvery}{imputation interval for saving imputed
  frequencies for the complete-data table. If \code{0}, then no
  imputations are saved.} 

  \item{saveProbSeries}{if \code{TRUE} then the simulated values of
  cell probabilities from MCMC will be stored as a series.}

  \item{typeMCMC}{either \code{"DA"} (data augmentation) or
  \code{"RWM"} (random-walk Metropolis); see DETAILS.}

  \item{tuneDA}{tuning parameters for data augmentation MCMC; see DETAILS.}

  \item{tuneRWM}{tuning parameter for random-walk Metropolis MCMC; see
  DETAILS.}

  \item{stuckLimit}{criterion for deciding if the MCMC algorithm has
  gotten stuck.}

  \item{startValDefault}{method used to obtain default starting values
  for parameters if no starting values are provided. \code{"center"}
  begins in the center of the parameter space, which assigns equal
  probability to all non-structural zero cells in the complete-data
  table. \code{"uniform"} draws
  random starting values from a uniform distribution on the cell
  probabilities.} 

  \item{startValJitter}{standard deviation for Gaussian random noise added to
  starting values. If \code{cvam} is called with
  \code{saturated=FALSE}, the log-linear coefficients are perturbed by
  this amount; if \code{saturated=TRUE}, the log-cell
  probabilities are perturbed by this amount and renormalized to sum
  to one.}

  \item{critEM}{convergence criterion for EM stopping rule; see DETAILS.}

  \item{critNR}{convergence criterion for Newton-Raphson stopping rule
  in M-step of EM; see DETAILS.}

  \item{critBoundary}{criterion for testing whether estimated probabilities
  are close to zero, in which case a warning is given.}

  \item{ncolMaxMM}{limit on the number of columns allowed for a
  log-linear model matrix.}

  \item{excludeAllNA}{if \code{TRUE}, then cases for which all modeled
  variables are missing will be excluded from the model fitting
  procedure, because they only contribute constant terms to the
  observed-data loglikelihood function.}

  \item{critModelCheck}{criterion for checking the log-linear model
  matrix for linear dependencies among the columns.}

  \item{confidence}{confidence coefficient for interval estimates,
  used when estimates are requested in the call to \code{cvam}.}

  \item{probRound}{if TRUE, estimated probabilities will be rounded.}

  \item{probDigits}{number of digits for rounding estimated probabilities.}

}

\details{
When \code{cvam} is called with \code{method="EM"}, it performs an EM
algorithm. At each E-step, observations with missing or coarsened
values are apportioned to cells of the complete-data table in the
expected amounts determined by the current estimated parameters. At
the M-step, the a log-linear model is fit to the predicted
complete-data frequencies from the E-step, using a Newton-Raphson procedure if
\code{saturated=FALSE}. The EM algorithm is stopped after
\code{iterMaxEM} iterations, or when the maximum
absolute difference in cell means from one iteration to the next is
no greater than \code{critEM}. The Newton-Raphson procedure in each
M-step is stopped after \code{iterMaxNR} iterations or when the
maximum absolute difference in cell means from one iteration to the next is
no greater than \code{critNR}.

When \code{cvam} is called with \code{method="MCMC"}, the algorithm
that is run depends on \code{typeMCMC} and on whether the model is fit
with \code{saturated=TRUE}.
\itemize{
\item If \code{saturated=FALSE} and
\code{typeMCMC="DA"}, then the algorithm is a data-augmentation
procedure that resembles EM. At each cycle, observations with missing
or coarsened values are randomly allocated to cells of the
complete-data table by drawing from a multinomial distribution, and
the log-linear coefficients are updated using one step of a
Metropolis-Hastings algorithm that mimics Newton-Raphson and
conditions on the allocated 
frequencies. The proposal distribution is multivariate-t and can be
adjusted by tuning constants in 
\code{tuneDA}, a numeric vector containing the degrees of
freedom, step size and scale factor.
\item If \code{saturated=FALSE} and \code{typeMCMC="RWM"},
the observations with missing or coarsened values are not allocated,
and the log-linear coefficients are updated by a step of random-walk
Metropolis. The proposal is mutivariate-t and can be adjusted by
tuning constants in \code{tuneRWM}, a numeric vector containing the
degrees of freedom and scale factor.
\item If \code{saturated=TRUE}, then the algorithm is a
data-augmentation procedure that requires no tuning.
}
Full details on the EM and MCMC procedures are given in the Appendix
of the vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.
}

\value{a list of control parameters for internal use by the
function \code{cvam}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{cvam}
}

\examples{
# display all control parameters and their default values
cvamControl()
}

