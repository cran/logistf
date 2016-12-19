\name{logistpl.control}
\Rdversion{1.1}
\alias{logistpl.control}
\title{Control Parameters for logistf Profile Likelihood Confidence Interval Estimation}
\description{
Sets parameters for modified Newton-Raphson iteration for finding profile likelihood confidence intervals
in Firth's penalized likelihood logistic regression
}
\usage{
logistpl.control(maxit=100, maxhs=5, maxstep=5, lconv=0.00001, xconv=0.00001, 
   ortho=FALSE, pr=FALSE)
}

\arguments{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a decrease in log likelihood.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood.}
  \item{xconv}{specifies the convergence criterion for the parameter estimates.}
  \item{ortho}{requests orthogonalization of variable for which confidence intervals are computed with respect to other covariates.}
  \item{pr}{request rotation of the matrix spanned by the covariates}
}
\details{
\code{logistpl.control()} is used by \code{logistf} to set control parameters to default values
when computing profile likelihood confidence intervals.
Different values can be specified, e. g., by \code{logistf(...,} 
\code{control=} \code{logistf.control(maxstep=1))}.
}


\value{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a decrease in log likelihood.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood.}
  \item{xconv}{specifies the convergence criterion for the parameter estimates.}
  \item{ortho}{specifies if orthogonalization is requested.}
  \item{pr}{specifies if rotation is requested}
}
\author{Georg Heinze}
\examples{
data(sexagg)
fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT, 
   plcontrol=logistpl.control(maxstep=1))
summary(fit2)
}
\keyword{regression}
\keyword{models}
