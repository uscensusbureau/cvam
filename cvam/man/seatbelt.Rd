\name{seatbelt}

\docType{data}

\alias{seatbelt}

\title{Seatbelt Data}

\description{
This dataset, previously analyzed by Hochberg (1977), Chen (1989) and
Schafer (1997), pertains to the effectiveness of seatbelts for 
preventing injury in automobile accidents. In a sample 
drawn from police reports, 80,084 accidents were classified by sex of
the driver, degree of damage to the car (low or high), seatbelt use
(no, yes) and whether the driver was injured (no, yes). Experience had
shown that police tended to overestimate the proportion of drivers who
were not using seatbelts and were not injured. To evaluate those
potential biases, a followup study was conducted with 1,796 additional
accidents. In the followup, investigators obtained more reliable
information on seatbelt use and injury from personal interviews and
hospital records.
}

\usage{seatbelt}

\format{a data frame with 80 rows and 8 variables:

\describe{    

    \item{\code{source}}{source of the data; \code{"sample"} indicates
    the original sample of police reports, and \code{"followup"}
    indicates the followup study}

    \item{\code{sex}}{sex of the driver (\code{"M"} or \code{"F"})}

    \item{\code{damage}}{degree of damage to the automobile
    (\code{"low"} or \code{"high"})}

    \item{\code{belt.p}}{whether the driver was wearing a seatbelt,
    according to the police report (\code{"no"} or \code{"yes"})}

    \item{\code{injury.p}}{whether the driver was injured,
    according to the police report (\code{"no"} or \code{"yes"})}

    \item{\code{belt.f}}{whether the driver was wearing a seatbelt,
    as determined in the followup study (\code{"no"} or \code{"yes"})}

    \item{\code{injury.f}}{whetherthe driver was injured,
    as determined in the followup study (\code{"no"} or \code{"yes"})}
}
}

\source{
Chen, T.T. (1989) A review of methods for misclassified categorical
data in epidemiology. \emph{Statistics in Medicine}, 8, 1095-1106.

Hochberg, Y. (1977) On the use of double sampling schemes in analyzing
categorical data with misclassification errors. \emph{Journal of the
American Statistical Association}, 72, 914-921.

Schafer, J.L. (1997) \emph{Analysis of Incomplete Multivariate
Data}. London: Chapman & Hall/CRC Press.
}

\keyword{ datasets }
