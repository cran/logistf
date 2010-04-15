\name{logistpl.control}
\alias{logistpl.control}
\title{Control parameters for logistf}
\description{
Sets parameters for Newton-Raphson iteration in Firth's penalized-likelihood logistic regression
}
\usage{
logistpl.control(maxit=100, maxhs=5, maxstep=5, lconv=0.00001, xconv=0.00001, ortho=FALSE, pr=FALSE)
}

\arguments{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a lower log likelihood than the target value.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood. For convergence, the loglikelihood may not be greater than target value + \code{lconv}}
  \item{xconv}{specifies the convergence criterion for the parameter estimates. For convergence, the change in any parameter estimate may not be greater than \code{xconv}}
  \item{ortho}{orthogonalizes the variable for which confidence intervals should be computed with respect to the other variables
               in the model. This may sometimes speed up or enable convergence. Default=FALSE}
  \item{pr}{replaces the matrix of 'other' variables in the model (those for which the current call of \code{logistpl} does not compute a
            confidence limit) by their principal components (of the same rank). This may sometimes speed up or enable convergence. Default=FALSE}
}
\value{
  \item{maxit}{the maximum number of iterations}
  \item{maxhs}{the maximum number of step-halvings in one iteration. The increment of the beta vector within one iteration is
               divided by 2 if the new beta leads to a decrease in log likelihood.}
  \item{maxstep}{specifies the maximum step size in the beta vector within one iteration.}
  \item{lconv}{specifies the convergence criterion for the log likelihood.}
  \item{xconv}{specifies the convergence criterion for the parameter estimates.}
  \item{ortho}{orthogonalizes the variable for which confidence intervals should be computed with respect to the other variables
               in the model. This may sometimes speed up or enable convergence. Default=FALSE}
  \item{pr}{replaces the matrix of 'other' variables in the model (those for which the current call of \code{logistpl} does not compute a
            confidence limit) by their principal components (of the same rank). This may sometimes speed up or enable convergence. Default=FALSE}
}
\details{
\code{logistpl.control()} is used by \code{logistf} (if \code{pl==TRUE}) and \code{logistfplot}, which both call the internal function \code{logistpl},
 to set control parameters to default values.
Different values can be specified, e. g., by \code{logistf(...,plcontrol=logistpl.control(maxstep=1})).}

\examples{
data(sexagg)
fit1<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT) # use default
fit1$pl.iter
fit2<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sexagg, weights=COUNT, plcontrol=logistpl.control(maxstep=1)) #change maxstep
fit2$pl.iter
}
\keyword{regression}
\keyword{models}
