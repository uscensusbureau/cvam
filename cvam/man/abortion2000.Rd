\name{abortion2000}

\docType{data}

\alias{abortion2000}

\title{Abortion Attitudes from the 2000 General Social Survey}

\description{This dataset, which was extracted from the 2000 General
Social Survey (GSS) (Smith et al., 2019), reports the responses of
adults in the United States to seven questions about legalized
abortion.  The questions began, ``Please tell me whether or not you
think it should be possible for a pregnant woman to obtain a legal
abortion if...'' The abortion items were given to a random two-thirds
subsample of GSS participants, so about 33\% of the values are \code{NA}
by design. Refusal to answer the question (a rare occurrence) was
also coded here as \code{NA}.

The data frame also includes variables on age, sex, race, Hispanic
origin, education, religious and political affiliation.  A second race
item, which was modeled after the race question on the U.S. Census
questionnaire, was given to a random half-sample.
}

\usage{abortion2000}

\format{a data frame with 2,817 rows and 16 factors:

\describe{    

    \item{\code{Age}}{respondent's age, with levels \code{"18-29"},
    \code{"30-49"}, \code{"50-64"}, and \code{"65+"}}

    \item{\code{Sex}}{respondent's sex, coded as \code{"Female"} or
    \code{"Male"}} 

    \item{\code{Race}}{respondent's race, coded as \code{"White"},
    \code{"Black"}, or \code{"Other"}; see NOTE below}

    \item{\code{CenRace}}{respondent's race, coded as \code{"White"},
    \code{"Black"}, \code{"Hisp"} or \code{"Other"}; see NOTE below}

    \item{\code{Hisp}}{respondent's Hispanic classification, with levels
    \code{"nonHisp"} and \code{"Hisp"}}

    \item{\code{Degree}}{respondent's education, classified as
    \code{"<HS"} (did not finish high school),
    \code{"HS"} (high school diploma), \code{"JunCol"} (junior
    college), \code{"Bach"} (Bacheor's degree), or \code{"Grad"}
    (graduate degree)} 

    \item{\code{Relig}}{respondent's religious preference, classified as
    \code{"Prot"} (Protestant),
    \code{"Cath"} (Roman Catholic), \code{"Jewish"}, \code{"None"}, or
    \code{"Other"}} 

    \item{\code{Party}}{respondent's political party identification, with levels
    \code{"Dem"} (Democrat), \code{"Rep"} (Republican), and
    \code{"Ind/Oth"} (Independent or Other); see NOTE below}

    \item{\code{PolViews}}{respondent's political views, with levels
    \code{"Con"} (Conservative), \code{"Mod"} (Moderate), and \code{"Lib"}
    (Liberal)}
}

Each of the seven variables below is a factor with levels
\code{"Yes"}, \code{"No"}, and \code{"DK"} (don't know). The items
 were prefixed by, ``Please tell me whether or not you think it
should be possible for a pregnant woman to obtain a legal abortion
if...''

\describe{

    \item{\code{AbDefect}}{``...If there is a strong chance of serious
    defect in the baby?''}

    \item{\code{AbNoMore}}{``...If she is married and does not want any
    more children?''}

    \item{\code{AbHealth}}{``...If the woman's own health is seriously
    endangered by the pregnancy?''}

    \item{\code{AbPoor}}{``...If the family has a very low income and
    cannot afford any more children?''}

    \item{\code{AbRape}}{``...If she became pregnant as a result of
    rape?''}

    \item{\code{AbSingle}}{``...If she is not married and does not
    want to marry the man''}

    \item{\code{AbAny}}{``...The woman wants it for any reason?''}

}
}

\note{

\code{Race}, which corresponds to the GSS variable \code{race}, is based
on the interviewer's assessment of the respondent's race. When
interviewers were not sure, they could ask the respondent, ``What race
do you consider yourself?''

\code{CenRace} is a collapsed version of the GSS variable
\code{racecen1}. That variable, which was modeled on the race question
in the U.S. Census, was given to half of the GSS sample in 2000 and
to the full sample in subsequent years. Participants could choose
from over a dozen race categories, or they could select ``Some other
race'' and provide their own. The \code{"Hisp"} values represent those
who chose ``Some other race'' and described themselves as Hispanic,
Latino, Latina, or something similar.

\code{Party} is based on the GSS variable \code{partyid}. Level
\code{"Dem"} includes Decocrat-leaning Independents, and 
\code{"Rep"} includes Republican-leaning Independents.

}

\source{
Smith, T.W., Davern, M., Freese, J., and Morgan, S.L. (2019)
\emph{General Social Surveys, 1972--2018}. National Data Program for
the Social Sciences, No. 25., 1 data file (64,814 logical records) + 1
codebook (3,758 pp.). Chicago: NORC.
}

\keyword{ datasets }
