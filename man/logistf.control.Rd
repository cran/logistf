\name{logistf.control}
\Rdversion{1.1}
\alias{logistf.control}
\title{Control parameters for logistf}
\description{
Sets parameters for Newton-Raphson iteration in Firth's penalized-likelihood logistic regression
}
\usage{
logistf.control(maxit = 25, maxhs = 5, maxstep = 5, lconv = 1e-05, gconv = 1e-05, 
    xconv = 1e-05, collapse=TRUE)
}

\arguments{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a decrease in log likelihood.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood.}
  \item{gconv}{specifies the convergence criterion for the first derivative of the log likelihood (the score vector).}
  \item{xconv}{specifies the convergence criterion for the parameter estimates.}
  \item{collapse}{if TRUE, evaluates all unique combinations of x and y and collapses data set. This may save computing time with large data sets with only categorical (binary) covariates.}
}
\details{
\code{logistf.control()} is used by \code{logistf} and \code{logistftest} to set control parameters to default values.
Different values can be specified, e. g., by \code{logistf(..., } \code{control=} \code{logistf.control(maxstep=1))}.}


\value{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a decrease in log likelihood.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood.}
  \item{gconv}{specifies the convergence criterion for the first derivative of the log likelihood (the score vector).}
  \item{xconv}{specifies the convergence criterion for the parameter estimates.}
  \item{collapse}{if TRUE, evaluates all unique combinations of x and y and collapses data set.}
}
\author{Georg Heinze}
\examples{
data(sexagg)
fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT, 
   control=logistf.control(maxstep=1))
summary(fit2)
}
\keyword{regression}
\keyword{models}
