\name{vcov.logistf}
\Rdversion{1.1}
\alias{vcov.logistf}
\title{
\code{vcov} Method for \code{logistf} Objects
}
\description{
Provides the \code{vcov} method for \code{logistf} objects.
}
\usage{
\method{vcov}{logistf}(object, ...)
}
\arguments{
  \item{object}{a \code{logistf}  object}
  \item{...}{further arguments}
}
\value{
the variance-covariance matrix of a \code{logistf} fit
}
\author{Georg Heinze}
\seealso{
\code{\link{logistf}}%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(sex2)
fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)
vcov(fit)
}
