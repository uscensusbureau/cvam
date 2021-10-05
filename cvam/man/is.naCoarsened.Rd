\name{is.naCoarsened}
\alias{is.naCoarsened}

\title{Missing-Value Indicators for Coarsened Factors}

\description{
A coarsened factor, produced by the function \code{\link{coarsened}},
stores missing values differently from an ordinary factor. If the
base R function \code{is.na} is applied to a coarsened factor, every element
of the result will be \code{FALSE}. The function \code{is.naCoarsened}
is the suitable alternative to \code{is.na} for coarsened factors.
}

\usage{
is.naCoarsened(x)
}

\arguments{
  \item{x}{a coarsened factor}
}

\details{
A coarsened factor, produced by the function \code{\link{coarsened}},
has two types of levels. Base levels represent states of complete
information, and coarse levels represent states of incomplete or
missing information. Each coarse level maps onto two or more base
levels. The last coarse level is \code{NA}, which
maps onto every base level.
} 

\value{
A logical vector of the same length as \code{x}, with \code{TRUE}
indicating that an element of \code{x} is \code{NA}, and \code{FALSE}
otherwise. 
}

\references{
For more information about coarsened factors in the \code{cvam}
package, see the vignette \emph{Understanding Coarsened Factors in cvam}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\seealso{
\code{\link{cvam}},
\code{\link{coarsened}},
\code{\link{baseLevels}},
\code{\link{dropCoarseLevels}}
}

\examples{
fac <- factor( c("red", "green", NA, "yellow", "notRed", "green") )
cFac <- coarsened( fac,
   levelsList = list("notRed" = c("green", "yellow")) )
is.naCoarsened(cFac)
}
