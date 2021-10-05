\name{cvam-package}
\alias{cvam-package}
\docType{package}
\title{
\packageTitle{cvam}
}
\description{
\packageDescription{cvam}
}
\details{

Log-linear models, when applied to frequencies in a multiway table,
describe associations among the factors that
define the dimensions of the table. Standard functions for fitting
log-linear models in R, including \code{\link{glm}} and
\code{\link{loglin}}, cannot accept observations with incomplete information
on any of the factors in the model, because those observations cannot
be assigned with certainty to a single cell of the complete-data table. The
functions in the \code{cvam} package facilitate log-linear modeling
of factors with missing and coarsened values. The two major functions
are: 
\preformatted{
   coarsened     Create a coarsened factor
   cvam          Log-linear models for coarsened factors
}
A coarsened factor is an extended version of a factor whose elements
may be fully observed, partially observed or missing.  The 
partially-observed and missing states are represented by extra levels
which are interpreted as groupings of the fully observed states.
The \code{cvam} function fits log-linear models to coarsened
factors. It also accepts ordinary factors with or without missing
values, and factors that contain only missing values,
which are useful for latent-class analysis. The modeling routines
implemented in \code{cvam} include EM algorithms for
mode finding and Markov chain Monte Carlo procedures for
Bayesian simulation and multiple imputation. Supporting
funtions are used to extract information from a fitted model,
including:
\preformatted{
   anova          Compare the fit of cvam models
   cvamEstimate   Estimated marginal and conditional probabilities
   cvamPredict    Predict missing or coarsened values
   cvamLik        Likelihood of observed data patterns
   cvamImpute     Impute missing or coarsened values
}
Five datasets are also provided:
\preformatted{
   abortion2000        Abortion attitudes from the General Social Survey
   crime               Crime victimization data
   hivtest             HIV test dataset
   microUCBAdmissions  U.C. Berkeley graduate admissions microdata
   seatbelt            Seatbelt data
}
}

\author{
\packageAuthor{cvam}

Maintainer: \packageMaintainer{cvam}
}

\references{
Extended descriptions and examples for all major functions are
provided in two vignettes, \emph{Understanding Coarsened Factors in
cvam} and \emph{Log-Linear Modeling with Missing and Coarsened Values
Using the cvam Package}.
}

\keyword{ package }
