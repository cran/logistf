\name{logistftest}
\alias{logistftest}
\alias{print.logistftest}
\title{Bias-reduced logistic regression}
\description{
This function performs a penalized likelihood ratio test on some (or
all) selected factors.  The resulting object is of the class logistftest
and includes the information printed by the proper print method.
}
\usage{
logistftest(formula=attr(data, "formula"), data=sys.parent(),
  test, values, maxit = 25, maxhs=5, epsilon = .0001,
  maxstep = 10, firth=TRUE, beta0)
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
  \item{test}{righthand formula of parameters to test (e.g. \code{~ B +
    D - 1}). As default all parameter apart from the intercept are tested.
    If -1 is not included in the formula, the intercept would be tested,
    too!  As alternative to the formula one can give the indexes of the
    ordered effects to test (a vector of integers). To test only the
    intercept specify \code{test = ~ - .} or \code{test = 1}.
  }
  \item{values}{null hypothesis values, default values are 0. For
    testing the specific hypothesis 1 = 1,4 = 2,5 = 0 we specify test= ~
    B1 + B4 + B5 - 1 and values=c(1, 2, 0).}
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
  \item{beta0}  {specifies the initial values of the coefficients for the fitting algorithm.}
}
\values{ 
The object returned is of the class logistf and has the following attributes:
  \item{testcov}{a vector of the fixed values of each covariate; NA stands for a parameter
    which is not tested.}
  \item{loglik}{a vector of the (penalized) log-likelihood of the full and the
    restricted models. If the argument beta0 not missing, the full model isn't
    evaluated.}
  \item{df: the number of degrees of freedom in the model.}
  \item{prob}{the p-value of the test.}
  \item{call}{the call object}
  \item{method}{depending on the fitting method "Penalized ML" or "Standard ML".}
  \item{beta}{the coefficients on the restricted solution.}

}
\details{ 
This function performs a penalized likelihood ratio test on some (or
all) selected factors.  The resulting object is of the class logistftest
and includes the information printed by the proper print method.

}

\references{
Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika} 
  80, 27--38.

Heinze G (1999). Technical Report 10: The application of Firth's procedure to Cox
  and logistic regression. Department of Medical Computer Sciences, Section of Clinical
  Biometrics, Vienna University, Vienna.

Heinze G, Schemper M (2002). A solution to the problem of 
  separation in logistic regression. \emph{Statistics in Medicine} 21: 2409-2419.

Heinze G, Ploner M (2003). Fixing the nonconvergence bug in 
logistic regression with SPLUS and SAS. \emph{Computer Methods and 
Programs in Biomedicine} 71: 181-187.

Ploner, M. (2001). Technical Report 2/2001: An SPLUS library to perform
  logistic regression without convergence problems. Section of Clinical Biometrics, Department of
  Medical Computer Sciences, University of Vienna, Vienna.
}
\seealso{logistf, logistfplot}

\examples{
data(sex2)
logistftest(case ~ age+oc+vic+vicl+vis+dia,  sex2, 
            test = ~ vic + vicl - 1, values = c(2, 0))
}

\keyword{regression}
\keyword{models}
