\name{summary.logistf}
\Rdversion{1.1}
\alias{summary.logistf}
\title{
\code{summary} Method for \code{logistf} Objects
}
\description{
Provides a summary method for \code{logistf}  objects.
}
\usage{
\method{summary}{logistf}(object, ...)
}
\arguments{
  \item{object}{a logistf object}
  \item{\dots}{additional parameters}
}
\value{
Side effect: a summary table with results from logistf analysis.
}
\author{Georg Heinze}
\examples{
data(sex2)
fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)
summary(fit)
}
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
