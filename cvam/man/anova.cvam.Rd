\name{anova.cvam}
\alias{anova.cvam}

\title{Comparing the Fit of Two or More Models}

\description{
Compares the fit of two or more \code{\link{cvam}} objects
}

\usage{
\method{anova}{cvam}(object, ..., method = c("lrt", "logP", "AIC", "BIC"),
   pval = FALSE, pvalDigits = 4L, showRank=NULL )
}

\arguments{
  \item{object}{an object produced by \code{cvam}}

  \item{\dots}{additional \code{cvam} objects}

  \item{method}{criterion for model comparison: \code{"lrt"} uses -2
  times the loglikelihood function for a standard likelihood-ratio test;
  \code{"logP"} uses -2 times the penalized loglikelihood or
  log-posterior density; \code{"AIC"} uses Akaike information criterion;
  and \code{"BIC"} uses Bayesian information criterion.}

  \item{pval}{if \code{TRUE} then p-values will be computed if method
  is \code{"lrt"} or \code{"logP"}.}

  \item{pvalDigits}{digits for rounding of p-values}

  \item{showRank}{if \code{TRUE}, models will be ranked from best to
  worst (with rank \code{1} being the best) according to the fit measure
  specified by \code{method}. Defaults to \code{TRUE} if \code{method} is
  \code{"AIC"} or \code{"BIC"}}
}

\details{
The p-values reported for the \code{"lrt"} and \code{"logP"} methods use
  a standard chi-squared approximation, with degrees of freedom equal
  to the difference in the number of parameters for the models being
  compared. This approximation is valid only if the models being
  compared are properly nested and ordered, with the simplest model
  appearing first in the argument list. The chi-squared approximation
  can be poor if the degrees of freedom for the comparison is large,
  and if the model corresponding to the null hypothesis (i.e., the
  smaller one) has fitted cell means that are too small. The
  chi-squared approximation is not appropriate for comparing
  latent-class models with a different number of latent classes.

The likelihood function used in \code{"lrt"} and \code{"logP"} is
based on a Poisson model for the cell means in the
complete-data table. The Poisson model is an appropriate surrogate for
a multinomial model. It is also an appropriate surrogate for a product
multinomial if the model includes all possible associations among
variables regarded as fixed.

The residual degrees of freedom are the difference between the number
of free parameters in a saturated Poisson model minus the number of free
parameters in the current Poisson model. The saturated model estimates one
parameter for every cell in the complete-data table, excluding
dimensions for latent factors, and excluding structural-zero
cells. No adjustments are made for estimates on a boundary of the
parameter space.

For \code{"BIC"}, the sample size is taken to be the total number of
observations or total frequency in the data supplied by the user to
fit the model, which does not include a flattening constant or any 
nuggets from a prior distribution created by
\code{\link{cvamPrior}}. No adjustments are made for missing or
coarsened values.
}

\value{an object of class \code{c("anova","data.frame")}}

\references{
For more information, refer to the package vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}}
}

\examples{
M0 <- cvam( ~ V1 + V2, data=crime, freq=n )
M1 <- cvam( ~ V1 * V2, data=crime, freq=n )
anova(M0, M1, pval=TRUE)
}
