\name{cig2019}

\docType{data}

\alias{cig2019}

\title{Cigarette Use from the 2019 National Health Interview Survey}

\description{The National Health Interview Survey (NHIS), conducted by
the National Center for Health Statistics and the U.S. Census Bureau,
is a large, annual, cross-sectional survey designed to measure health
status and behaviors in the non-institutionalized
U.S. population. This dataset, which was extracted from the 2019 NHIS,
includes responses to several questions pertaining to 
cigarettes. Adult participants (age 18+) were asked, ``Have
you smoked at least 100 cigarettes in your entire life?'' Those
who answered ``Yes'' were asked additional questions about
current use; for those who did not, the interviewers skipped to the
next section.

Analyses of NHIS data should account for the stratified multistage
cluster sample design. Sample weights, stratum and cluster indicators
are included for this purpose.
}

\usage{cig2019}

\format{a data frame with 31,997 rows and 14 variables:

\describe{    

    \item{\code{wtia_a}}{interim annual weight for analysis of adult
    questionnaire items; see NOTE below}

    \item{\code{wtfa_a}}{final annual weight for analysis of adult
    questionnaire items; see NOTE below}

    \item{\code{urbrrl}}{urban/rural classification scheme for
    counties, a factor with levels \code{"Large central metro"},
    \code{"Large fringe metro"}, \code{"Med/small metro"}. and
    \code{"Nonmetropolitan"}}
    
    \item{\code{region}}{region within the United States, a factor
    with levels \code{"Northeast"}, \code{"Midwest"}, \code{"South"},
    and \code{"West"}}
    
    \item{\code{pstrat}}{pseudo-stratum for variance estimation}
    
    \item{\code{ppsu}}{pseudo-primary sampling unit for variance
    estimation; see NOTE below}
    
    \item{\code{agep_a}}{age in years; see NOTE below}
    
    \item{\code{sex_a}}{sex, a factor with levels \code{"Male"} and
    \code{"Female"}}
    
    \item{\code{hispall_p}}{single and multiple race groups with
    Hispanic origin, a factor with levels \code{"Hispanic"}, \code{"NH
    White only"}, \code{"NH Black/AfrAmer only"}, \code{"NH Asian
    only"}, \code{"NH AIAN only"}, \code{"NH AIAN and other group"},
    and \code{"Other single and multiple races"}}
    
    \item{\code{smkev_a}}{``Have you smoked at least 100 cigarettes in
    your ENTIRE LIFE?'', a factor with levels \code{"Yes"} and
    \code{"No"}} 
    
    \item{\code{smknow_a}}{``Do you NOW smoke cigarettes every day,
    some days or not at all?'', a factor with levels \code{"Every
    day"}, \code{"Some days"}, and \code{"Not at all"}; see NOTE below} 
    
    \item{\code{cignow_a}}{``On average, about how many cigarettes do
    you NOW smoke a day?''; see NOTE below}
    
    \item{\code{smk30d_a}}{``On how many of the PAST 30 DAYS did you
    smoke a cigarette?''; see NOTE below}
    
    \item{\code{cig30d_a}}{``On average, when you smoked during the
    PAST 30 DAYS, about how many cigarettes did you smoke a day?'';
    see NOTE below}
}
}

\note{
The weight variables \code{wtia_a} and \code{wtfa_a} both reflect the
probability that the person was selected into the sample. The
interim annual adult weight \code{wtia} includes adjustments for 
nonresponse due to non-participation in the interview. The final
annual adult weight \code{wtia} includes adjustments for nonresponse
and additional factors to calibrate the sample to population totals
for important geographic and demographic variables.

The cluster identifier \code{ppsu} is nested within the stratum
identifier \code{pstrat}; for example, \code{ppsu == 1} in
\code{pstrat == 100} and \code{ppsu == 1} in \code{pstrat == 101}
represent different clusters.

The age variable \code{agep_a} is top-coded at 85; the value
\code{agep_a == 85} indicates 85 or older.

Question \code{smknow_a}, ``Do you NOW smoke cigarettes every day,
some days, or not at all?'', was asked only if \code{smkev_a =="Yes"};
otherwise it is missing.

Question \code{cignow_a}, ``On average, about how many cigarettes
do you NOW smoke a day?'', was asked only if \code{smknow_a == "Every
day"}; otherwise it is missing. This variable is top-coded at \code{95}.

Question \code{smk30d_a}, ``On how many of the PAST 30 DAYS did
you smoke a cigarette?'', was asked only if \code{smknow_a == "Some
days"}; otherwise it is missing.

Question \code{cig30d_a}, ``On average, when you smoked during the
PAST 30 DAYS, about how many cigarettes did you smoke a day?'', was
asked only if \code{smknow_a == "Some days"} and if
\code{(smk30d_a > 0) | is.na(smk30d_a)}. This variable is top-coded at
\code{95}. 

}

\source{
National Center for Health Statistics. National Health Interview
Survey, 2019. Public-use data file and documentation. 
\href{https://www.cdc.gov/nchs/nhis/data-questionnaires-documentation.htm}{https://www.cdc.gov/nchs/nhis/data-questionnaires-documentation.htm}. 2020.
}

\keyword{ datasets }
