\name{add1.logistf}
\alias{add1.logistf}
\alias{drop1.logistf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add or Drop All Possible Single Terms to a logistf Model
}
\description{
Compute all the single terms in the scope argument that can be added to or dropped from the model, fit those models and compute a table of the changes in fit.
}
\usage{
\method{add1}{logistf}(object, scope, test = "PLR", ...)
\method{drop1}{logistf}(object, scope, test = "PLR", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted logistf object}
  \item{scope}{
The scope of variables considered for adding or dropping. Should be a vector of variable names. Can be left missing; the method will then use all variables in the object`s data slot which are not identified as the response variable.}
  \item{test}{
The type of test statistic. Currently, only the PLR test (penalized likelihood ratio test) is allowed for logistf fits.
}
  \item{\dots}{
further arguments passed to or from other methods.}
}
\details{
\code{drop1} and \code{add1} generate a table where for each variable the penalized likelihood ratio chi-squared, the degrees of freedom, and the p-value for dropping/adding this variable are given.
}
\value{
A matrix with \code{nvar} rows and 3 columns (Chisquared, degrees of freedom, p-value).
}
\author{
Georg Heinze
}
\examples{
data(sex2)
fit<-logistf(data=sex2, case~1, pl=FALSE)
add1(fit)

fit2<-logistf(data=sex2, case~age+oc+dia+vic+vicl+vis)
drop1(fit2)
}
\keyword{ models }
\keyword{ regression }
