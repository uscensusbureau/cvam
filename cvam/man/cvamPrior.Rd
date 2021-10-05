\name{cvamPrior}
\alias{cvamPrior}
\alias{summary.cvamPrior}
\alias{summary.cvamPrior}
\alias{print.cvamPrior}
\alias{print.summary.cvamPrior}

\title{Data-Augmentation Prior for Coarsened Factor Loglinear Model}

\description{
The \code{\link{cvam}} function fits loglinear models to coarsened
categorical variables. The \code{cvamPrior} function creates an object
to pass to \code{cvam} to represent prior information that is incorporated
into the model fit.}

\usage{
cvamPrior(obj = list(), flatten = 0, ridge = 0, intensity = 1)
}

\arguments{
  \item{obj}{a list of prior information nuggets to apply to the
  complete-data frequency table; see DETAILS.}
  \item{flatten}{a prior information nugget to be divided equally
  across all cells of the complete-data frequency table; see DETAILS.}
  \item{ridge}{a ridge factor to apply to the log-linear coefficients;
  see DETAILS.}
  \item{intensity}{a factor applied simultaneously to all prior information
  to scale it up or down; see DETAILS.}
}

\details{
An object produced by this function, when passed to \code{cvam}
through its \code{prior} argument, incorporates prior information as
\itemize{
\item a flattening constant, a positive value that is divided equally
among all non-structural zero cells of the complete-data table, and
\item prior nuggets, which take the form of coarsened-data frequencies
that are assigned to selected cells or groups of cells.
}
Log-linear models fit with \code{saturated=FALSE} can also accept a
ridge factor, 
which acts upon the coefficients in a manner similar
to ridge regression, shrinking the estimated coefficients 
toward zero and stabilizing its estimated
covariance matrix. The added information is equivalent to
a multivariate normal prior density centered at zero
with prior precision (inverse covariance) matrix equal to the ridge
factor times the identity matrix.

The intensity factor provides a simple way to strengthen or weaken
the overall amount of prior information, which is useful for
sensitivity analyses. The flattening constant, nugget frequencies and
ridge factor are all multiplied by \code{intensity}. Setting
\code{intensity=2} doubles the prior information, \code{intensity=.5}
cuts it in half, and so on.
}

\value{an object of class \code{"cvamPrior"}, designed for use by the
function \code{\link{cvam}}.
}

\references{
For more information, refer to the package vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}},
\code{\link{coarsened}}
}

\examples{
# fit a saturated model to a four-way table
fit <- cvam( ~ Sex*CenRace*Hisp*Party, data=abortion2000,
 saturated=TRUE )
# add a flattening constant
fit <- cvam( ~ Sex*CenRace*Hisp*Party, data=abortion2000,
 saturated=TRUE, prior=cvamPrior( flatten=10 ) )

# fit with saturated=FALSE and no prior information, and
# notice how large the SEs are
fit <- cvam( ~ Sex*CenRace*Hisp*Party, data=abortion2000,
 saturated=FALSE )
head( get.coef(fit, withSE=TRUE) )
# add a very mild ridge factor and notice how the SEs
# have become reasonable
fit <- cvam( ~ Sex*CenRace*Hisp*Party, data=abortion2000,
 saturated=FALSE, prior=cvamPrior( ridge=.1 ) )
head( get.coef(fit, withSE=TRUE) )

# add s few prior nuggets to stabilize the distribution
# of Party within a rare category
nuggetList <- list(
   list( CenRace="Black", Hisp="Hisp", Party="Dem", freq=1 ),
   list( CenRace="Black", Hisp="Hisp", Party="Rep", freq=1 ),
   list( CenRace="Black", Hisp="Hisp", Party="Ind/Oth", freq=1 ) )
myPrior <- cvamPrior( nuggetList, flatten=10 )
summary(myPrior)
fit <- cvam( ~ Sex*CenRace*Hisp*Party, data=abortion2000,
 saturated=FALSE, prior=myPrior )
}
