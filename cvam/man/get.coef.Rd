\name{get.coef}
\alias{get.coef}
\alias{get.covMat}
\alias{get.estimates}
\alias{get.loglik}
\alias{get.logP}
\alias{get.mfTrue}
\alias{get.modelMatrix}
\alias{get.offset}
\alias{get.strZero}
\alias{get.fitted}
\alias{get.imputedFreq}
\alias{get.minus2logPSeries}
\alias{get.coefSeries}
\alias{get.probSeries}

\title{Extract Information from a Coarsened-Variable Model}

\description{This group of functions will extract various
summaries from a model fit by \code{\link{cvam}}, either
with the EM algorithm or by Markov chain Monte Carlo.
}

\usage{

get.coef(obj, withSE = FALSE, meanSeries = TRUE, msgIfNone = TRUE)

get.covMat(obj, msgIfNone = TRUE)

get.estimates(obj, msgIfNone = TRUE)

get.loglik(obj, msgIfNone = TRUE)

get.logP(obj, msgIfNone = TRUE)

get.mfTrue(obj)

get.modelMatrix(obj, msgIfNone = TRUE)

get.offset(obj, mfTrue = FALSE, msgIfNone = TRUE)

get.strZero(obj, mfTrue = FALSE, msgIfNone = TRUE)

get.fitted(obj, type=c("prob", "mean", "logMean"), mfTrue = TRUE, 
   meanSeries = TRUE, msgIfNone = TRUE )

get.imputedFreq(obj, msgIfNone = TRUE)

get.minus2logPSeries(obj, startValShift = TRUE,
   msgIfNone = TRUE, coda = ( obj$method == "MCMC" ) )

get.coefSeries(obj, msgIfNone = TRUE, coda = ( obj$method == "MCMC" ) )

get.probSeries(obj, levelNames=TRUE, sep=".",
msgIfNone = TRUE, coda = ( obj$method == "MCMC" ) )

}


\arguments{

\item{obj}{an object resulting from a call to \code{\link{cvam}}
  with \code{method = "EM"} or \code{method = "MCMC"}}

\item{withSE}{if \code{TRUE}, then \code{get.coef} will return a data
frame containing estimated log-linear coefficients, standard errors,
t-statistics and p-values; if \code{FALSE}, then only a vector of
coefficients is given.}

\item{mfTrue}{if \code{TRUE}, then \code{get.offset},
\code{get.strZero} and \code{get.fitted} will return a data frame with
one row per cell, with all model variables (the non-coarsened
versions) present as factors, and with another variable named
(depending on which function was called) \code{offset}, \code{strZero}
or \code{fit} containing the requested values. If \code{FALSE}, then 
\code{get.offset}, \code{get.strZero} and \code{get.fitted} will
return a vector containing the requested values.}

\item{meanSeries}{applies when \code{obj} is the result from a simulation
run. If \code{TRUE}, results will be based on from a running average
of simulated parameters over all iterations after the burn-in
period. If \code{FALSE}, results will be based only on the simulated
parameter values at the end of the run.}

\item{msgIfNone}{if \code{TRUE} then, if the \code{get} procedure
fails, an informative message is given 
explaining why the requested summary cannot be obtained. For example,
\code{get.coef} will fail to return coefficients from a model fit with
\code{cvam(..., saturated = TRUE)} because no model matrix is created
and the log-linear coefficients
are not defined. If \code{FALSE}, then these messages are suppressed.}

\item{type}{type of fitted values to be returned by
\code{get.fitted}. \code{"prob"} returns cell probabilities
conditioned on variables fixed by the model (if any); \code{"mean"}
returns cell means from the log-linear model; and \code{"logMean"}
returns log-cell means from the log-linear model (i.e., the linear
predictor).}

\item{startValShift}{the function \code{get.minus2logPSeries} extracts a
saved series from an MCMC run containing  the values of (minus 2 times)
the log-posterior density function. If \code{startValShift} is true,
the series is shifted by (minus 2 times) the log-posterior density at
the starting value, if the starting value appears to be a mode.}

\item{coda}{if \code{TRUE}, the series from an MCMC run is returned as
an \code{mcmc} object for plotting and diagnostic analysis with the
\code{coda} package. If \code{FALSE}, a one-dimensional series is
returned as a numeric vector, and a multidimensional series is
returned as a numeric matrix with rows corresponding to iterations and
columns corresponding to elements of the multidimensional quantities
being monitored.}

\item{levelNames}{the \code{get.probSeries} function extracts a saved
series of probabilities from an MCMC run corresponding to cells of the
complete-data table (i.e., the rows of \code{mfTrue}). If
\code{levelNames} is \code{TRUE}, names for the cell probabilities are
constructed from the levels of the factors in \code{mfTrue}. As the
number of variables in the model grows, these names can become
unwieldly, and setting \code{levelNames} to \code{FALSE} omits the
names.}

\item{sep}{a character string used to separate the levels of multiple
factors when \code{levelNames} is \code{TRUE}.}

}

\details{

The series objects returned by \code{get.minus2logPSeries},
\code{get.coefSeries} and \code{get.probSeries} omit results from the
burn-in period, if any, and may also be thinned. The default behavior
is no burn-in period and no thinning. The burn-in period and thinning
interval are set by components of the \code{control} argument to
\code{cvam}, via the function \code{\link{cvamControl}}; the relevant
components are \code{control$burnMCMC} and \code{control$thinMCMC}.
By default, \code{cvam} does not save cell probabilities. To save them, set
\code{control$saveProbSeries} to \code{TRUE}.

\code{get.imputedFreq} returns multiple imputations of frequencies for
the complete-data table generated and stored during an MCMC run after
the burn-in period. The default behavior is no imputation. This can be
changed by setting \code{control$imputeEvery} to an integer greater
than zero.

Other useful information from a model fit can be extracted
with the \code{summary} method for a \code{cvam} object, and with the
functions \code{\link{cvamEstimate}}, \code{\link{cvamPredict}},
\code{\link{cvamLik}}, and \code{\link{cvamImpute}}.
}

\value{

\code{get.coef} returns a vector of estimated coefficients from the
log-linear model if \code{withSE=FALSE}; if \code{withSE=TRUE}, it
returns a data frame containing coefficients, standard errors,
t-statistics and p-values.

\code{get.covMat} returns an estimated covariance matrix for the
estimated coefficients.

\code{get.estimates} returns a data frame or a list of data frames
containing the estimates held in \code{obj}.

\code{get.loglik} and \code{get.logP} return the value of the
loglikelihood function or log-posterior density from the beginning of
the final iteration of EM or MCMC. If the model was fit using
\code{cvam(..., saturated=TRUE)}, the likelihood is based on a
multinomial or product-multinomial distribution over the cells of the
complete-data table. If the model was fit as a log-linear approach
using \code{cvam(..., saturated=FALSE)}, the likelihood is based on a
surrogate Poisson model.

\code{get.mfTrue} returns a data frame with one row per cell of the 
complete-data table. The variables in this data frame include every
factor appearing in the model (the non-coarsened versions) and another
variable named \code{freq}. If the model was fit using \code{cvam(...,
method="EM")}, \code{freq} contains the predicted cell
frequencies at the final iteration of EM. If the model was fit using
\code{cvam(..., method="MCMC")}, \code{freq} contains a running
average of imputed cell frequencies over all iterations of MCMC after the
burn-in period. In either case, if the data used to fit the model
contain no missing or coarsened values, then \code{freq} will be equal
to the observed frequencies.

\code{get.modelMatrix} returns the model matrix for the log-linear
model. The rows of the model matrix correspond to the rows of
\code{mfTrue}, and the columns correspond to terms created from
the factors in \code{mfTrue}.

\code{get.offset} retrieves the offset for the log-linear model. 
If \code{mfTrue} is \code{TRUE}, it returns the data frame
\code{mfTrue} with a numeric variable named \code{offset}. If
\code{mfTrue} is \code{FALSE}, it returns a numeric vector of length
\code{NROW(mfTrue)}.

\code{get.strZero} retrieves the logical values indicating whether
each cell is structural zero. If \code{mfTrue} is \code{TRUE}, it
returns the data frame \code{mfTrue} with a logical variable named
\code{strZero}. If \code{mfTrue} is \code{FALSE}, it returns a logical
vector of length \code{NROW(mfTrue)}.

\code{get.fitted} retrieves fitted values from the log-linear
model. If \code{type="prob"}, the fitted values are cell
probabilities conditioned on any variables fixed in the model. If
\code{type="mean"} or \code{"logMean"}, the fitted values are cell
means or log-cell means from the log-linear model.
If \code{mfTrue} is \code{TRUE}, the function
returns the data frame \code{mfTrue} with a numeric variable named
\code{fit}. If \code{mfTrue} is \code{FALSE}, it returns a numeric
vector of length \code{NROW(mfTrue)}.

\code{get.imputedFreq} returns the data frame \code{mfTrue}, with the
\code{freq} variable replaced by multiply imputed versions of the
frequencies for the complete-data table.

\code{get.minus2logPSeries} returns a series of (minus 2 times) the
log-posterior density values from the iterations of MCMC, either as a
numeric vector or as an \code{mcmc} object used by the \code{coda}
package. 

\code{get.coefSeries} returns a series of log-linear coefficients from
the iterations of MCMC, either as a numeric matrix or as an
\code{mcmc} object used by the \code{coda} package.

\code{get.probSeries} returns a series of cell probabilities from
the iterations of MCMC, either as a numeric matrix or as an
\code{mcmc} object used by the \code{coda} package.

}

\references{
For more information, refer to the package vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.

For information about \code{coda}, see:

Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). CODA:
Convergence Diagnosis and Output Analysis for MCMC, \emph{R News},
vol 6, 7-11.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}}, \code{\link{summary.cvam}}, \code{\link{cvamControl}},
\code{\link{cvamEstimate}}, \code{\link{cvamPredict}},
\code{\link{cvamImpute}}, \code{\link{cvamLik}}
}

\examples{
fit <- cvam( ~ V1 * V2, data=crime, freq=n )
get.coef(fit, withSE=TRUE)
get.covMat(fit)
get.fitted(fit, type="mean")

set.seed(6755)
fit <- cvam(fit, method="MCMC",
   control=list(iterMCMC=5000, imputeEvery=500) )
get.imputedFreq(fit)

\dontrun{plot( get.coefSeries(fit) )  # coda trace and density plots}
}
