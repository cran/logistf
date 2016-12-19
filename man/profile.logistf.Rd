\name{profile.logistf}
\alias{profile.logistf}
\alias{plot.profile.logistf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Profile Penalized Likelihood
}
\description{
Evaluates the profile penalized likelihood of a variable based on a logistf model fit.
}
\usage{
\method{profile}{logistf}(fitted, which, variable, steps = 100, pitch = 0.05, limits, 
   alpha = 0.05, firth = TRUE, legends = TRUE, control, plcontrol, plot = FALSE,
   ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitted}{
An object fitted by \code{logistf}
}
  \item{which}{
A righthand formula to specify the variable for which the profile should be evaluated, e.g., \code{which=~X}).
}
  \item{variable}{
Alternatively to \code{which}, a variable name can be given, e.g., \code{variable="X"}
}
  \item{steps}{
number of steps in evaluating the profile likelihood
}
  \item{pitch}{
alternatively to \code{steps}, one may specify the step width in multiples of standard errors
}
  \item{limits}{
lower and upper limits of parameter values at which profile likelihood is to be evaluated
}
  \item{alpha}{
the significance level (1-\eqn{\alpha} the confidence level,
    0.05 as default).
}
  \item{firth}{
use of Firth's penalized maximum likelihood (\code{firth=TRUE}, default) or the
    standard maximum likelihood method (\code{firth=FALSE}) for the logistic regression. 
}
  \item{legends}{
legends to be included in the optional plot
}
  \item{control}{
Controls Newton-Raphson iteration. Default is \code{control=} \code{logistf.control(maxstep,}
 \code{maxit, maxhs, lconv, gconv, xconv})
}
  \item{plcontrol}{
Controls Newton-Raphson iteration for the estimation of the profile likelihood confidence intervals.
       Default is \code{plcontrol=} \code{logistpl.control(maxstep,} 
       \code{maxit,}  \code{maxhs, lconv, xconv, ortho, pr})
}
  \item{plot}{
If TRUE, profile likelihood is plotted. This parameter becomes obsolete as a generic \code{plot} function is now provided.
}
  \item{\dots}{
Further arguments to be passed.}
}
\value{
An object of class \code{logistf.profile} with the following items:
\item{beta}{parameter values at which likelihood was evaluated}
\item{stdbeta}{parameter values divided by standard error}
\item{profile}{profile likelihood, standardized to 0 at maximum of likelihood. The values in profile are given as minus \eqn{\chi^2}.}
\item{loglike}{unstandardized profile likelihood}
\item{signed.root}{signed root (\eqn{z}) of \eqn{\chi^2} values (negative for values below the maximum likelihood estimate, positive for values above the maximum likelihood estimate)}
\item{cdf}{profile likelihood expressed as cumulative distribution function, obtained as \eqn{\Phi(z)}, where \eqn{\Phi} denotes the standard normal distribution function.}
}
\author{
Georg Heinze and Meinhard Ploner}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{plot.profile.logistf}}
 }
\examples{
data(sex2)
fit<-logistf(case ~ age+oc+vic+vicl+vis+dia, data=sex2)
plot(profile(fit,variable="dia"))
plot(profile(fit,variable="dia"), "cdf")
plot(profile(fit,variable="dia"), "density")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
