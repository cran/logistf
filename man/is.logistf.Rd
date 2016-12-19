\name{is.logistf}
\alias{is.logistf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check \code{logistf} Objects
}
\description{
Checks if an object is of class \code{logistf}
}
\usage{
is.logistf(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object.
}
}
\value{
A logical value
}
\author{
Georg Heinze
}
\examples{
data(sex2)
fit<-logistf(data=sex2, case~1, pl=FALSE)
is.logistf(fit)

}
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
