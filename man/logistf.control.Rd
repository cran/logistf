\name{logistf.control}
\alias{logistf.control}
\title{Control parameters for logistf}
\description{
Sets parameters for Newton-Raphson iteration in Firth's penalized-likelihood logistic regression
}
\usage{
logistf.control(maxit=25, maxhs=5, maxstep=5, lconv=0.00001, gconv=0.00001, xconv=0.00001)
}

\arguments{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a decrease in log likelihood.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood.}
  \item{gconv}{specifies the convergence criterion for the first derivative of the log likelihood (the score vector).}
  \item{xconv}{specifies the convergence criterion for the parameter estimates.}
}
\value{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a decrease in log likelihood.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood.}
  \item{gconv}{specifies the convergence criterion for the first derivative of the log likelihood (the score vector).}
  \item{xconv}{specifies the convergence criterion for the parameter estimates.}
}
\details{
\code{logistf.control()} is used by \code{logistf} and \code{logistftest} to set control parameters to default values.
Different values can be specified, e. g., by \code{logistf(...,control=logistf.control(maxstep=1})).}

\examples{
data(sexagg)
fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT, control=logistf.control(maxstep=1))
summary(fit2)
}
\keyword{regression}
\keyword{models}
