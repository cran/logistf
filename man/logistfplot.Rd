\name{logistfplot}
\alias{logistfplot}
\title{Plot penalized profile likelihood}
\description{
This function plots the penalized profile likelihood for a specified
parameter.
}
\usage{
logistfplot <- function(formula = attr(data, "formula"),
    data = sys.parent(), which, pitch = 0.05, limits, alpha = 0.05,
    maxit = 25, maxhs = 5, epsilon = 0.0001, maxstep = 10, firth = TRUE, legends = TRUE)
}

\arguments{
  \item{formula}{a formula object, with the response on the left of the  operator, and the
    model terms on the right. The response must be a vector with 0 and 1 or FALSE and
    TRUE for the model outcome, where the higher value (1 or TRUE) is modeled. It's possible
    to include contrasts, interactions, nested effects, cubic or polynomial splines and all
    the S-PLUS features, as well, e.g. \code{Y ~ X1^*X2 + ns(X3, df=4)}.
  }
  \item{data}{a data.frame where the variables named in the formula can be found, i. e.
    the variables containing the binary response and the covariates.
  }
  \item{which}{a righthand formula specifying the plotted parameter, interaction or
    general term, e.g. \code{~ A - 1} or \code{~ A : C - 1}. The profile likelihood of the 
    intercept would be obtained by the formula \code{~ - .}.}
  \item{pitch}{distances between the interpolated points in standard errors of
    the parameter estimate, the default value is 0.05.}
  \item{limits}{vector of the minimum and the maximum on the x-scale in standard
    deviations distant form the maximum likelihood. The default values
    are the extremes of both confidence intervals, Wald and PL, plus or minus half a
    standard deviation of the parameter, respectively.}
  \item{alpha}{the significance level (1-\eqn{\alpha} the confidence level,
    0.05 as default).}
  \item{maxit}{maximum number of iterations (default value is 25)}
  \item{maxhs}{maximum number of step-halvings per iterations (default value is 5)}
  \item{epsilon}{specifies the maximum allowed change in penalized log likelihood to
    declare convergence. Default value is 0.0001.}
  \item{maxstep}{specifies the maximum change of (standardized) parameter values allowed
    in one iteration. Default value is 0.5.}
  \item{firth}{use of Firth's penalized maximum likelihood (firth=TRUE, default) or the
    standard maximum likelihood method (firth=FALSE) for the logistic regression. Note
    that by specifying pl=TRUE and firth=FALSE (and probably a lower number of iterations) 
  one obtains profile likelihood confidence intervals for maximum likelihood logistic
  regression parameters.}
  \item{\beta0}  {specifies the initial values of the coefficients for the fitting algorithm.}
  \item{legends}{if FALSE, legends on the bottom of the plot would be omitted
    (default is TRUE).}
}
\values{ 
The object returned is a simple data.frame containing three columns which
allow reproducing the plot. Each row represents one point of the interpolation. The
columns are as follows:
  \item{std}{distance from the maximum of the profile likelihood (in standard
   errors of the parameter estimate).}
  \item{name}{the value of the parameter for the variable name specified
   in argument \code{which}.}
  \item{loglik.pen}{the value of the penalized likelihood.}
}
\details{ 
This function plots the profile likelihood of a specific parameter based
on the penalized likelihood.  A symmetric shape of the profile penalized
log likelihood (PPL) function allows use of Wald intervals, while an
asymmetric shape demands profile penalized likelihood intervals (Heinze
& Schemper (2001)).


}

\references{
Heinze G (1999). Technical Report 10: The application of Firth's procedure to Cox
  and logistic regression. Department of Medical Computer Sciences, Section of Clinical
  Biometrics, Vienna University, Vienna.

Heinze G, Schemper M (2001). A solution to the problem of monotone likelihood
  in logistic regression. submitted.

}
\seealso{logistf, logistftest}

\keyword{regression}
\keyword{models}
