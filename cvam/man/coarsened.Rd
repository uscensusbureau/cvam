\name{coarsened}
\alias{coarsened}
\alias{is.coarsened}
\alias{print.coarsened}
\alias{droplevels.coarsened}
\alias{relevel.coarsened}
\alias{reorder.coarsened}
\alias{rep.coarsened}
\alias{[.coarsened}
\alias{[[.coarsened}
\alias{[<-.coarsened}
\alias{[[<-.coarsened}

\title{Coarsened Factors}

\description{A coarsened factor is 
an extended version of a factor or ordered factor whose elements
may be fully observed, partially observed or missing.  The 
partially-observed and missing states are represented by extra levels
which are interpreted as groupings of the fully observed states.
Coarsened factors are specifically designed for modeling with the
\code{cvam} package.
}

\usage{
coarsened(obj, levelsList = list(), warnIfCoarsened = TRUE)

is.coarsened(x)

\method{print}{coarsened}(x, quote = FALSE, max.levels = NULL,
   width = getOption("width"), ...)

\method{droplevels}{coarsened}(x, ...)

\method{relevel}{coarsened}(x, ...)

\method{reorder}{coarsened}(x, ...)

\method{rep}{coarsened}(x, ...)

\method{[}{coarsened}(x, ...)

\method{[[}{coarsened}(x, ...)

\method{[}{coarsened}(x, ...) <- value

\method{[[}{coarsened}(x, ...) <- value
}

\arguments{
  \item{obj}{a factor or ordered factor to be converted to a coarsened
  factor}

  \item{levelsList}{a named list that defines the groupings of
  \code{levels(obj)} to indicate states of partial knowledge}

  \item{warnIfCoarsened}{if \code{TRUE}, a warning is issued if
  \code{obj} is already a coarsened factor}

  \item{x}{a coarsened factor or other object}

  \item{quote}{logical, indicating whether or not strings should be
  printed with surrounding quotes}

  \item{max.levels}{integer, indicating how many base levels and
  coarse levels should be printed for a coarsened factor; if \code{0}, no
  extra base levels or coarse levels lines will be printed. The
  default, \code{NULL}, entails choosing \code{max.levels} such that 
  the base levels and coarse levels each print on one line of width
  \code{width}}

  \item{width}{only used when \code{max.levels} is \code{NULL}; see
  above}

  \item{...}{additional arguments passed to or from other methods}

  \item{value}{character: a set of levels for replacement}
}

\details{
A coarsened factor, which inherits from class \code{"factor"} or
\code{c("ordered", "factor")}, has two types of levels: base
levels, which represent states of complete knowledge, and
coarse levels, which represent states of incomplete
knowledge. Each coarse level maps to two or more
base levels. The mapping is defined by the argument \code{levelsList}.

For example, consider a factor whose levels are 
\code{c("red", "notRed", "green", "yellow")}, where \code{"notRed"} denotes an
observation that is either \code{"green"} or \code{"yellow"}. When the
factor is converted to a coarsened factor, \code{c("red", "green",
"yellow")}  becomes the \code{\link{baseLevels}}, and \code{"notRed"} becomes
an element of \code{\link{coarseLevels}}. To produce this result, the
argument \code{levelsList} should have a component named
\code{"notRed"}, whose value is \code{c("green", "yellow")}.

The last coarse level is \code{NA}, denoting an observation
that could belong to any of the base levels. The \code{NA} coarse level is
created automatically. Calling \code{coarsened} with an empty
\code{levelsList} (the default) produces a coarsened factor with
\code{NA} as its only coarse level. 

If the main argument to \code{coarsened} is already a coarsened
factor, then a warning is issued (if \code{warnIfCoarsened} is
\code{TRUE}) and the coarsened factor is returned unchanged.

The generic functions \code{\link{droplevels}}, \code{\link{relevel}},
and \code{\link{reorder}} should not be applied to coarsened factors;
the S3 methods \code{droplevels.coarsened}, \code{relevel.coarsened},
and \code{reorder.coarsened} will prevent this from happening.

\code{rep.coarsened} is a method for the generic function \code{\link{rep}}
that ensures the special attributes of a coarsened factor are preserved.

Extraction and replacement methods \code{`[`} and \code{`[[`} are also
provided to preserve the special attributes of coarsened factors.
}

\value{
\code{coarsened} returns a coarsened factor.

\code{is.coarsened} returns \code{TRUE} if \code{x} is a coarsened
factor and \code{FALSE} otherwise.
}

\references{
For more information, refer to the package vignette
\emph{Understanding Coarsened Factors in cvam}.
}

\author{Joe Schafer \email{Joseph.L.Schafer@census.gov} }

\note{
Coarsened factors were designed for use by the modeling function
\code{\link{cvam}}, which treats base levels and coarse levels
differently. Other statistical modeling routines, such as \code{lm},
may not handle them appropriately. Functions outside of the
\code{cvam} package will treat coarse levels (including \code{NA}) the
same as base levels, producing results that are difficult to
interpret or nonsensical, especially if the base levels are ordered.

The behavior of \code{coarsened} with \code{levelsList = list()} is
similar to that of \code{\link{addNA}}, which converts the missing
values in a factor to non-missing observations with value \code{NA}
and adds \code{NA} to the levels. The result of \code{addNA}, however,
is an ordinary factor or ordered factor which has no mechanism to
inform other functions that \code{NA} has special meaning.

The function \code{is.na} should not be applied to a coarsened factor;
use \code{\link{is.naCoarsened}} instead.

Because base levels and coarse levels should be handled differently,
functions from base R that manipulate the levels of a factor,
including \code{relevel}, \code{reorder}, \code{droplevels}, and
the replacement version of \code{levels} should not
be used with coarsened factors. Supplying a coarsened factor to any of
these functions will produce an error.
}

\seealso{
\code{\link{cvam}},
\code{\link{is.naCoarsened}},
\code{\link{baseLevels}},
\code{\link{dropCoarseLevels}}
}

\examples{
fac <- factor( c("red", "green", NA, "yellow", "notRed", "green") )
cFac <- coarsened( fac,
   levelsList = list("notRed" = c("green", "yellow")) )
print(cFac)
# extraction and replacement
print( cFac[2:3] )
cFac[2:3] <- c("NA", "green") 
}
