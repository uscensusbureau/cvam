\name{cvam}
\alias{cvam}
\alias{cvam.formula}
\alias{cvam.cvam}
\alias{cvam.default}

\title{
Log-Linear Models for Incomplete Categorical Variables 
}

\description{
Fits log-linear models to categorical variables by three
methods: maximizing the loglikelihood or log-posterior density by
Expectation-Maximization (EM) algorithms, simulating the posterior
distribution by a Markov chain Monte Carlo (MCMC) algorithms, and
creating random draws of parameters from an approximate Bayesian
posterior distribution. The factors in the model may have missing
or coarsened values.
}

\usage{

% the generic function
cvam(obj, ...)

\method{cvam}{formula}(obj, data, freq, prior = cvamPrior(),
    method = c("EM", "MCMC", "approxBayes", "mfSeen", "mfTrue",
       "mfPrior", "modelMatrix"), control = list(), omitData = FALSE,
    saturated = FALSE, modelMatrix = NULL, offset = NULL,
    strZero = NULL, startVal = NULL, estimate = NULL, ...)

\method{cvam}{cvam}(obj, method = obj$method, control = NULL, startVal = NULL, 
    estimate = NULL, ...)

}

\arguments{
  \item{obj}{an object used to select a method: either a model
  formula or the result from a previous call to \code{cvam}.}

  \item{data}{an optional data frame, list or environment (or object
  coercible to a data frame by \code{as.data.frame}) containing the variables
  in the model. If not found in \code{data}, the variables are taken from
  \code{environment(obj)}, typically the environment from which
  \code{cvam} is called.}

  \item{freq}{an optional variable for holding integer frequencies when the
  observations are grouped. If \code{freq} is not given, then the
  observations are assumed to represent microdata, and all frequencies
  are set to one.}

  \item{prior}{an object produced by \code{\link{cvamPrior}} to
  represent prior information incorporated into the model fit.}

  \item{method}{a procedure for fitting the model:
  \code{"EM"} computes a maximum-likelihood (ML) estimate, penalized
  ML estimate, or posterior mode; \code{"MCMC"} runs a Markov chain
  Monte Carlo algorithm to simulate a sequence of correlated random
  draws from the posterior distribution of the unknown parameters;
  \code{"approxBayes"} creates independent draws from an approximate
  posterior distribution. The other alternatives return various
  objects without fitting the model.}

  \item{control}{a named list containing control parameters which are
  passed to \code{\link{cvamControl}}. Control parameters determine the
  maximum number of iterations, criteria for judging convergence,
  proposal distributions for MCMC, and so on. Control parameters that are
  not found in this list are set to default values.}

  \item{omitData}{if \code{TRUE}, then the observations supplied
  through \code{data} and \code{freq} are ignored, and the fitted model
  is based only the prior information supplied through
  \code{prior}. Combining \code{omitData=TRUE} with
  \code{method="MCMC"} will simulate random draws from the prior
  distribution.}

  \item{saturated}{if \code{TRUE}, then a saturated model is fit to
  the cell means without defining a model matrix or log-linear
  coefficients.}
  
  \item{modelMatrix}{an optional model matrix that defines the
  log-linear model. In ordinary circumstances, \code{cvam} creates the
  model matrix automatically by interpreting terms in the model
  formula and referring to the \code{contrast} attributes of the
  model factors. In rare circumstances, a user may want to supply a
  different model matrix. The model matrix should have one row for
  every cell in the complete-data table. If a model matrix is
  supplied, the model formula is used only to identify the variables
  that are included the model, not to define the associations among them.}

  \item{offset}{an optional numeric vector of length
  \code{NROW(modelMatrix)} containing an offset for the log-linear
  model. If omitted, the offset is assumed to be zero for every cell.}
  
  \item{strZero}{an optional logical vector of length
  \code{NROW(modelMatrix)} containing
  \code{TRUE} for every cell to be considered a structural zero and
  \code{FALSE} elsewhere. Structural zeros are assumed to have zero
  probability and are omitted from the model fitting. If
  \code{strZero} is omitted, all elements are assumed to be \code{FALSE}.}

  \item{startVal}{an optional vector of starting values for the model
  parameters. If \code{saturated=FALSE}, this should be a vector of
  length \code{NCOL(modelMatrix)} containing log-linear coefficients;
  if \code{saturate=FALSE}, it should be a vector of length
  \code{NROW(modelMatrix)} containing cell probabilities or cell
  means, which are automatically rescaled to become probabilities.}

  \item{estimate}{an optional formula or list of formulas of the kind
  expected by \code{\link{cvamEstimate}} specifying marginal or
  conditional probabilities to be estimated, bypassing the need for a
  subsequent call to that function.}

  \item{\dots}{values to be passed to the methods.}
}

\details{
A log-linear model is specified by a one-sided formula that determines
which associations among the variables are allowed. For
example, \code{~ A + B + C} implies that \code{A}, \code{B} and
\code{C} are mutually independent; \code{~ A*B + A*C} implies that
\code{B} and \code{C} are conditionally independent given \code{A};
and so on. Variables in a model may be factors or coarsened factors,
and missing values are permitted. All models are fit using a surrogate
Poisson formulation which is appropriate for Poisson, multinomial
or product-multinomial sampling. A formula may contain a vertical bar
to specify variables to be regarded as fixed; for example, \code{~ A*B
+ A*C | A} fixes the variable \code{A}. Fixing variables does not
change the model fitting procedure; the only difference is that, after
the model has been fit, the cell probabilities are scaled to sum to
one within every combination of levels of the fixed variables.

If \code{cvam} is called with a \code{cvam} object as its first
argument, then the data, model and prior distribution will be
taken from the previous run, and (unless \code{startVal} is
supplied), starting values will be set to the final parameter values
from the previous run. 
}

\value{if \code{method} is \code{"EM"}, \code{"MCMC"} or
\code{"approxBayes"}, an object of class \code{c("cvam","list")}
containing the results of a model fit. For other values of
\code{method}, the requested object is returned without fitting a
model.}

\references{
Extended descriptions and examples for all major functions are
provided in two vignettes, \emph{Understanding Coarsened Factors in
cvam} and \emph{Log-Linear Modeling with Missing and Coarsened Values
Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{coarsened}},
\code{\link{cvamPrior}},
\code{\link{cvamControl}},
\code{\link{cvamEstimate}},
\code{\link{get.coef}},
\code{\link{summary.cvam}}
}


\examples{
# convert U.C. Berkeley admissions three-way table to data frame,
# fit model of conditional independence, display summary
# compare the fit to the saturated model
dF <- as.data.frame(UCBAdmissions)
fit <- cvam( ~ Dept*Gender + Dept*Admit, data=dF, freq=Freq )
summary(fit)
fitSat <- cvam( ~ Dept*Gender*Admit, data=dF, freq=Freq )
anova(fit, fitSat, pval=TRUE)

# fit non-independence model to crime data; then run MCMC for
# 5000 iterations, creating 10 multiple imputations of the frequencies
# for the 2x2 complete-data table
fit <- cvam( ~ V1 * V2, data=crime, freq=n )
set.seed(56182)
fitMCMC <- cvam(fit, method="MCMC", 
   control=list( iterMCMC=5000, imputeEvery=500) )
get.imputedFreq(fitMCMC)
}
