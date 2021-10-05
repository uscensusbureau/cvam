\name{baseLevels}
\alias{baseLevels}
\alias{nBaseLevels}
\alias{coarseLevels}
\alias{nCoarseLevels}
\alias{mapping}

\title{Get Coarsened Factor Attributes}

\description{Retrieve specific attributes of a coarsened factor.}

\usage{
baseLevels(x)
nBaseLevels(x)
coarseLevels(x)
nCoarseLevels(x)
mapping(x)
}

\arguments{
  \item{x}{a coarsened factor}
}

\details{
A coarsened factor, produced by the function \code{\link{coarsened}},
is an extended type of factor whose elements may be fully observed,
partially observed, or missing. The full set of attributes of a
coarsened factor may be obtained by \code{attributes}, and
individual attributes are available with \code{attr}. The functions
documented on this page are convenient alternatives to \code{attr}.

The elements of \code{baseLevels}, a character vector of length
\code{nBaseLevels}, represent states of complete knowledge. The
elements of \code{coarseLevels}, a character vector of length
\code{nCoarseLevels}, represent states of incomplete or missing
information. Each element of \code{coarseLevels} maps onto two or more
elements of \code{baseLevels}. The attribute \code{mapping} is an
integer matrix with \code{nCoarseLevels} rows and \code{nBaseLevels}
columns, with \code{1} in position \code{[i,j]} if coarse level
\code{i} contains base level \code{j}. The last coarse level is
always \code{NA}, and it contains every base level.
} 

\value{
The requested attribute of \code{x}.
}

\references{
For more information, refer to the package vignette
\emph{Understanding Coarsened Factors in cvam}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\note{
A coarsened factor has the usual attributes of a factor, but they
should not be altered directly.  For example, 
the function \code{levels<-}, the replacement version
of \code{levels}, should not be used with a coarsened factor.
}

\seealso{
\code{\link{cvam}},
\code{\link{coarsened}},
\code{\link{is.naCoarsened}},
\code{\link{dropCoarseLevels}}
}

\examples{
fac <- factor( c("red", "green", NA, "yellow", "notRed", "green") )
cFac <- coarsened( fac,
   levelsList = list("notRed" = c("green", "yellow")) )
baseLevels(cFac)
mapping(cFac)
}
