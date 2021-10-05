\name{cvamLik}
\alias{cvamLik}

\title{
Likelihood of Observed Data Patterns
}
\description{
After fitting a log-linear model with \code{cvam}, 
the fitted model object may be passed to this function, along with a
dataset that may contain missing or coarsened values, to compute
the likelihood of each pattern of possibly incomplete or
coarsened data for subset of variables,
possibly conditioned upon another subset of variables
}

\usage{
cvamLik(form, obj, data, meanSeries = TRUE)
}

\arguments{
  \item{form}{a formula indicating which variables to consider, and
  which variables to condition on, when computing the likelihood}

  \item{obj}{an object produced by \code{\link{cvam}} containing 
  results from a model fit}

  \item{data}{data frame for computing the likelihood values,
  possibly different from the data used to fit the model contained in
  \code{obj}}

  \item{meanSeries}{applies when \code{obj} contains results from a
     simulation run.
     If \code{TRUE}, then the requested likelihood values are based on
     a running mean of cell probabilities over all iterations 
     after the burn-in period. If \code{FALSE}, then the requested
     values are based only on the cell probabilities from the final
     iteration, and (assuming the run was sufficiently long, if MCMC)
     can be regarded as a single draw from their posterior distribution.} 

}

\value{
A data frame containing the model variables, with a variable
\code{likVall} holding the likelihood values
}

\details{
For structural zeros, \code{0/0} is returned as \code{0}. If any
variables are being conditioned on in \code{form}, they must not
contain any missing or coarsened values.
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
\code{\link{cvamPredict}}
}

\examples{
result <- cvam( ~ V1 * V2, freq=n, data=crime)
cvamLik( ~ V1 + V2, result, data=crime )
}
