\name{hivtest}

\docType{data}

\alias{hivtest}

\title{HIV test dataset}

\description{
This dataset concerns the diagnostic accuracy of tests for HIV
infection.  Four different tests (\code{A}, \code{B},
\code{C}, and \code{D}) were applied to 428 high-risk patients.  The
result from each test is either negative (\code{"neg"}) or positive
(\code{"pos"}). Yang and Becker (1997) applied latent-class analysis
to these data to estimate the sensitivity and specificity of these
tests in the absence of a gold standard.
}

\usage{hivtest}

\format{a data frame with 9 rows and 5 variables:

\describe{    

    \item{\code{A}}{test result from radioimmunoassay (RIA) utilizing
    recombinant agl21}

    \item{\code{B}}{test result from RIA utilizing purified HIV p24}

    \item{\code{C}}{test result from RIA utilizing purified HIV
    gpl20}

    \item{\code{D}}{test result from enzyme-linked immunosorbent assay} 

    \item{\code{COUNT}}{frequency; number of patients exhibiting the given 
    pattern of results}

}
}

\source{
Yang, I and Becker, M.P. (1997) Latent variable modeling of diagnostic
accuracy. \emph{Biometrics}, 53, 948-958.
}

\keyword{ datasets }
