\name{dropCoarseLevels}
\alias{dropCoarseLevels}

\title{Remove Coarse Levels from a Coarsened Factor}

\description{
A coarsened factor, produced by the function \code{\link{coarsened}},
has two types of levels: base levels, to represent observations that
are fully known, and coarse levels, to represent observations that are
partially or completely missing. The function \code{dropCoarseLevels}
converts a coarsened factor to a factor or ordered factor by
removing all of the coarse levels and setting the
corresponding observations to \code{NA}.
}

\usage{
dropCoarseLevels(x)
}

\arguments{
  \item{x}{a factor or coarsened factor}
}

\details{
If the only coarse level of \code{x} is \code{NA}, then no information
is lost when \code{dropCoarseLevels} is applied. If \code{x} has
other non-empty coarse levels, then the partial information carried by
those observations is effectively discarded.
} 

\value{
A factor or ordered factor, obtained by removing the coarse levels of
\code{x} and setting the observations in those levels to \code{NA}.

If \code{x} is a factor but not a coarsened factor, then it is
returned unchanged. 
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
\code{\link{is.naCoarsened}}
}

\examples{
fac <- factor( c("red", "green", NA, "yellow", "notRed", "green") )
cFac <- coarsened( fac,
   levelsList = list("notRed" = c("green", "yellow")) )
dropCoarseLevels(cFac)
}
