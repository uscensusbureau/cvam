\name{microUCBAdmissions}

\docType{data}

\alias{microUCBAdmissions}

\title{UC Berkeley Graduate Admissions Microdata}

\description{The dataset \code{UCBAdmissions}, distributed with the R
\code{datasets} package, is a three-dimensional array that 
classifies applicants to graduate school at the University
of California, Berkeley in 1973 by sex, department, and admission
status. That table, first published by Bickel \emph{et al.} (1975),
is frequently used in textbooks and statistics courses to illustrate
Simpson's paradox. This dataset, \code{microUCBAdmissions} is a
simulated data frame with three factors and one row per individual;
when the factors are tabulated, they reproduce the frequencies in
\code{UCBAdmissions}.
}

\usage{microUCBAdmissions}

\format{a data frame with 4,526 rows and 3 factors:

\describe{    

    \item{\code{Admit}}{\code{"Admitted"} or \code{"Rejected"}}

    \item{\code{Gender}}{\code{"Male"} or \code{"Female"}} 

    \item{\code{Dept}}{\code{"A"}, \code{"B"}, \code{"C"}, \code{"D"},
    \code{"E"}, or \code{"F"}}

}
}

\source{
Bickel, P.J., Hammel, E.A. and O'Connell, J.W. (1975) Sex bias in
graduate admissions: Data from Berkeley. \emph{Science}, 187, 398-403.
}

\keyword{ datasets }
