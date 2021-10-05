\name{summary.cvam}
\alias{print.cvam}
\alias{summary.cvam}
\alias{print.summary.cvam}

\title{Summarize a \code{cvam} object}

\description{Summarizes the result from a call to \code{cvam}}

\usage{
\method{summary}{cvam}(object, showCoef=TRUE, showEstimates=TRUE,
   digits = 4L, ...)
\method{print}{cvam}(x, ...)
\method{print}{summary.cvam}(x, ...)
}

\arguments{
  \item{object}{an object resulting from a call to \code{cvam}
  containing results  from a fitted model.}
  \item{showCoef}{if \code{TRUE}, the table of coefficients will be
  displayed when the \code{print} method is invoked.}
  \item{showEstimates}{if \code{TRUE}, estimated marginal and
  conditional probabilities requested by the \code{estimate} argument
  will be displayed when the \code{print} method is invoked.}
  \item{digits}{for printing.}
  \item{...}{additional arguments to be passed to methods.}
  \item{x}{\code{cvam} or \code{cvam.summary} object}
}

\value{for the \code{summary} method, an object of class
\code{"summary.cvam"}.}

\references{
For more information, refer to the package vignette \emph{Log-Linear
Modeling with Missing and Coarsened Values Using the cvam Package}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }


\seealso{
\code{\link{cvam}}
}

\examples{
# saturated model for hivtest data
result <- cvam( ~ A*B*C*D, data=hivtest, freq=COUNT )
summary(result)
}
