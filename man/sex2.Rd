\name{sex2}
     \docType{data}
     \alias{sexagg}
     \alias{sex2}
     \title{Urinary tract infection in American college students}
     \description{
       This data set deals with urinary tract infection in sexually active college women, along with covariate
       information on age an contraceptive use.
       The variables are all binary and coded in 1 (condition is present) and 0 (condition is absent): case (urinary tract infection, the study outcome variable),
       age (>= 24 years), dia (use of diaphragm), oc (use of oral contraceptive), vic (use of condom),
      vicl (use of lubricated condom), and vis (use of spermicide).
     }
     \usage{
     sex2
     sexagg}
     \format{\code{sex2}: a \code{data.frame} containing 239 observations; \code{sexagg}: an aggregated \code{data.frame} containing 31 observations with case weights (\code{COUNT}).}
     \source{www.cytel.com, the full data was available for download at their webpage, and is distributed along with their LogXact software. Here a reduced version of the original data is reproduced.}
     \references{
       Cytel Inc., (2010) \emph{LogXact 9 user manual}, Cambridge, MA:Cytel Inc.
     }
     \keyword{datasets}