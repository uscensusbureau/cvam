\name{latentFactor}
\alias{latentFactor}
\alias{is.latentFactor}

\title{Latent Factor}

\description{
A latent factor is a categorical variable whose values are entirely
missing. The function \code{latentFactor} provides a convenient way to
create a latent factor with a given number of base levels, which is
useful for latent-class modeling with \code{\link{cvam}}.
}

\usage{
latentFactor( n, levels = 2L )

is.latentFactor(x)
}

\arguments{
  \item{n}{length of the factor}

  \item{levels}{either an integer specifying the number of base levels, or
  a character vector containing labels for the base levels}

  \item{x}{an object to be tested}
}


\value{
For \code{latentFactor}, a latent coarsened factor of length
\code{n}; for \code{is.latentFactor}, \code{TRUE} or \code{FALSE}.
}

\references{
For more information, refer to the package vignettes
\emph{Understanding Coarsened Factors in  
cvam} and \emph{Log-Linear Modeling with Missing and Coarsened Values
Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}},
\code{\link{coarsened}},
\code{\link{baseLevels}},
\code{\link{is.naCoarsened}}
}

\examples{
# fit latent class model to hivtest data
hivtest$L <- latentFactor( NROW(hivtest), 2 )
set.seed(125)
fit <- cvam( ~ L*A + L*B + L*C + L*D, data=hivtest, freq=COUNT,
   control=list(startValJitter=.1) )
cvamEstimate( list( ~L, ~A|L, ~B|L, ~C|L, ~D|L ), fit )
}
