% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profile.logistf.r
\encoding{UTF-8}
\name{profile.logistf}
\alias{profile.logistf}
\title{Compute Profile Penalized Likelihood}
\usage{
\method{profile}{logistf}(
  fitted,
  which,
  variable,
  steps = 100,
  pitch = 0.05,
  limits,
  alpha = 0.05,
  firth = TRUE,
  legends = TRUE,
  control,
  plcontrol,
  ...
)
}
\arguments{
\item{fitted}{An object fitted by \code{logistf}}

\item{which}{A righthand formula to specify the variable for which the profile should be evaluated, e.g., which=~X).}

\item{variable}{Alternatively to which, a variable name can be given, e.g., variable="X"}

\item{steps}{Number of steps in evaluating the profile likelihood}

\item{pitch}{Alternatively to steps, one may specify the step width in multiples of standard errors}

\item{limits}{Lower and upper limits of parameter values at which profile likelihood is to be evaluated}

\item{alpha}{The significance level (1-\eqn{\alpha} the confidence level, 0.05 as default).}

\item{firth}{Use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default)
or the standard maximum likelihood method (\code{firth=FALSE}) for the logistic regression.}

\item{legends}{legends to be included in the optional plot}

\item{control}{Controls Newton-Raphson iteration. Default is \code{control= logistf.control(maxstep, 
maxit, maxhs, lconv, gconv, xconv)}}

\item{plcontrol}{Controls Newton-Raphson iteration for the estimation of the profile likelihood
confidence intervals. Default is \code{plcontrol= logistpl.control(maxstep, maxit, maxhs, lconv, xconv, ortho, pr)}}

\item{...}{Further arguments to be passed.}
}
\value{
An object of class \code{logistf.profile} with the following items:
\item{beta}{Parameter values at which likelihood was evaluated}
\item{stdbeta}{Parameter values divided by standard error}
\item{profile}{profile likelihood, standardized to 0 at maximum of likelihood. The values in
profile are given as minus \eqn{\chi^2}}
\item{loglik}{Unstandardized profile likelihood}
\item{signed.root}{signed root (z) of \eqn{\chi^2} values (negative for values below the maximum likelihood
estimate, positive for values above the maximum likelihood estimate)}
\item{cdf}{profile likelihood expressed as cumulative distribution function, obtained as
\eqn{\Phi(z)}, where \eqn{\Phi} denotes the standard normal distribution function.}
}
\description{
Evaluates the profile penalized likelihood of a variable based on a logistf model fit
}
\examples{
data(sex2)
fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)
plot(profile(fit,variable="dia"))
plot(profile(fit,variable="dia"), "cdf")
plot(profile(fit,variable="dia"), "density")

}
\references{
Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining
profile likelihood information from logistic regressions. Statistics in Medicine, to appear.
}
