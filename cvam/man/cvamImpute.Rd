\name{cvamImpute}
\alias{cvamImpute}


\title{Impute Data from a Fitted Model}

\description{
After fitting a log-linear model with \code{cvam}, 
the fitted model object may be passed to this function, along with a
dataset containing missing or coarsened values, to randomly impute the
true data from their predictive distribution given the
observed data and parameters from the fitted model.
}

\usage{
cvamImpute(obj, data, freq, meanSeries = FALSE, synthetic=FALSE)
}

\arguments{

  \item{obj}{an object produced by \code{\link{cvam}} containing 
  results from a model fit}

  \item{data}{data frame for imputation, possibly different from the
  data used to fit the model contained in \code{obj}}

  \item{freq}{variable containing frequencies for 
  \code{data}. If omitted, all frequencies are taken to be \code{1},
  meaning that the imputation frame is assumed to contain microdata.}

  \item{meanSeries}{applies when \code{obj} contains results from a
     simulation run. If \code{TRUE}, then the imputations are based on
     a running mean of cell probabilities over all iterations 
     after the burn-in period. If \code{FALSE}, then the imputations
     are based only on the cell probabilities from the final
     iteration, and (assuming the run was sufficiently long, if MCMC)
     can be regarded as a single draw from their posterior distribution. 
     See DETAILS below.}

  \item{synthetic}{if \code{TRUE}, then observed values for all
  variables in the data frame (excluding variables that are conditioned on
  in the model and regarded as fixed) are set to \code{NA} and 
  imputed, producing a dataset that is fully synthetic.}
 
}

\note{
When this function is used within a process for multiple imputation,
\code{meanSeries} should 
be set to \code{FALSE}, otherwise the imputations will not correctly reflect
uncertainty about model parameters. 
}

\value{a data frame containing imputed data. If \code{freq} was given,
the data frame has one row for each cell in the complete-data table
and a variable \code{freq} containing the frequencies. If \code{freq}
was not given, the data frame has one row for each microdata observatiom.
}

\references{
For more information, refer to the package vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}},
\code{\link{cvamEstimate}},
\code{\link{cvamPredict}},
\code{\link{cvamLik}}
}

\examples{
# impute from a grouped dataset with frequencies
fit <- cvam( ~ V1 * V2, freq=n, data=crime )
cvamImpute( fit, data=crime, freq=n )
# impute microdata
fit <- cvam( ~ Sex * PolViews * AbAny, data=abortion2000 )
impData <- cvamImpute( fit, data=abortion2000 )
head(impData)
}
