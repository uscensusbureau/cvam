\name{crime}

\docType{data}

\alias{crime}

\title{Crime Victimization Data}

\description{

This dataset, reported and analyzed by Kadane (1985), comes from the
National Crime Survey conducted by the U.S. Bureau of the
Census. Occupants of housing units were interviewed to determine
whether they had been victimized by crime in the preceding six-month
period. Six months later, the units were visited again to determine
whether the occupants had been victimized during the intervening
months. Missing values for various reasons occurred at both occasions.
}

\usage{crime}

\format{a data frame with 9 rows and 3 variables:

\describe{    

    \item{\code{V1}}{victimization status from the first visit}

    \item{\code{V2}}{victimization status from the second visit}

    \item{\code{n}}{frequency in sample}

}
}

\source{
Kadane, J.B. (1985) Is victimization chronic? a Bayesian analysis of
multinomial missing data. \emph{Journal of Econometrics}, 29, 47-67.

Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
}

\keyword{ datasets }
