\name{cvamPredict}
\alias{cvamPredict}

\title{Predict Missing or Coarsened Values from a Fitted Model}

\description{
After fitting a log-linear model with \code{cvam}, 
the fitted model object may be passed to this function, along with a
dataset containing missing or coarsened values, to predict one or more
variables from their predictive distribution given the
observed data and parameters from the fitted model.
}

\usage{
cvamPredict(form, obj, data, freq, meanSeries = TRUE, sep = ".")
}

\arguments{
  \item{form}{a one-sided formula indicating the variable or
  variables to be predicted, with variables separated by `\code{+}'}

  \item{obj}{an object produced by \code{\link{cvam}} containing 
  results from a model fit}

  \item{data}{data frame for prediction, possibly different from the
  data used to fit the model contained in \code{obj}}

  \item{freq}{variable containing frequencies for 
  \code{data}. If omitted, all frequencies are taken to be \code{1},
  meaning that the prediction frame is assumed to contain microdata.}

  \item{meanSeries}{applies when \code{obj} contains results from a
     simulation run.
     If \code{TRUE}, then the requested predictions are based on
     a running mean of cell probabilities over all iterations 
     after the burn-in period. If \code{FALSE}, then the requested
     predictions are based only on the cell probabilities from the final
     iteration, and (assuming the run was sufficiently long, if MCMC)
     can be regarded as a single draw from their posterior distribution.} 

  \item{sep}{character sting used to separate the levels of multiple
  variables being predicted}
}

\details{
Predictions from this function are unlike predictions from a
regression model. In regression, prediction is to compute the
estimated mean response ar specific values of the predictors. With
this function, predictions are based on the predictive distribution for
one or more variables given all the observed data, including the
variable(s) to be predicted if they are seen. The prediction for a
variable that is seen will assign a probability of one to the seen
value and zero probability to other values.
}

\value{
A data frame containing the predicted probabilities or frequencies,
with an attribute \code{colFrame} that identifies its columns
}

\references{
For more information, refer to the package vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}},
\code{\link{cvamEstimate}},
\code{\link{cvamImpute}},
\code{\link{cvamLik}}
}

\examples{
fit <- cvam( ~ V1 + V2, freq=n, data=crime )
cvamPredict( ~ V1, fit, data=crime, freq=n )   # predict frequencies
cvamPredict( ~ V1, fit, data=crime  )          # predict probabilities
}
