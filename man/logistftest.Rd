\name{logistftest}
\Rdversion{1.1}
\alias{logistftest}
\title{Penalized Likelihood Ratio Test}
\description{
This function performs a penalized likelihood ratio test on some (or
all) selected factors.  The resulting object is of the class logistftest
and includes the information printed by the proper print method.
}
\usage{
logistftest(object, test, values, firth = TRUE, beta0, weights, control)
}

\arguments{
  \item{object}{a fitted logistf object
  }
  \item{test}{righthand formula of parameters to test (e.g. \code{~ B +
    D - 1}). As default all parameter apart from the intercept are tested.
    If the formula includes -1, the intercept is omitted from testing.
    As alternative to the formula one can give the indexes of the
    ordered effects to test (a vector of integers). To test only the
    intercept specify \code{test = ~ - .} or \code{test = 1}.
  }
  \item{values}{null hypothesis values, default values are 0. For
    testing the specific hypothesis B1=1, B4=2, B5=0 we specify \code{test= } 
    \code{~B1+B4+B5-1} and \code{values=c(1,} \code{2,0)}.}
  \item{firth}{use of Firth's (1993) penalized maximum likelihood (\code{firth=TRUE}, default) or the
    standard maximum likelihood method (\code{firth=FALSE}) for the logistic regression. Note
    that by specifying \code{pl=TRUE} and \code{firth=FALSE} (and probably a lower number of iterations)
    one obtains profile likelihood confidence intervals for maximum likelihood logistic
    regression parameters.}
  \item{beta0}{specifies the initial values of the coefficients for the fitting algorithm.}
  \item{weights}{case weights}
  \item{control}{control parameters for iterative fitting}

}
\details{
This function performs a penalized likelihood ratio test on some (or
all) selected factors.  The resulting object is of the class logistftest
and includes the information printed by the proper print method.  Further
documentation can be found in Heinze & Ploner (2004).

In most cases, the functionality of the logistftest function is replaced by \code{anova.logistf}, which is a more 
standard way to perform likelihood ratio tests. However, as shown in the example below, \code{logistftest} provides
some specials such as testing agains non-zero values. (By the way, \code{anova.logistf} calls  \code{logistftest}.)

A \code{print} method is available.
}
\value{
The object returned is of the class logistf and has the following attributes:
  \item{testcov}{a vector of the fixed values of each covariate; NA stands for a parameter
    which is not tested.}
  \item{loglik}{a vector of the (penalized) log-likelihood of the full and the
    restricted models. If the argument beta0 not missing, the full model isn't
    evaluated.}
  \item{df}{the number of degrees of freedom in the model.}
  \item{prob}{the p-value of the test.}
  \item{call}{the call object}
  \item{method}{depending on the fitting method `Penalized ML' or `Standard ML'.}
  \item{beta}{the coefficients on the restricted solution.}
}
\references{
Firth D (1993). Bias reduction of maximum likelihood estimates. \emph{Biometrika}
  80, 27--38.

Heinze G, Ploner M (2004). Technical Report 2/2004: A SAS-macro, S-PLUS library and
  R package to perform logistic regression without convergence problems.
  Section of Clinical Biometrics, Department of
  Medical Computer Sciences, Medical University of Vienna, Vienna, Austria.
   \url{http://www.meduniwien.ac.at/user/georg.heinze/techreps/tr2_2004.pdf}

Heinze G (2006). A comparative investigation of methods for logistic regression
with separated or nearly separated data. \emph{Statistics in Medicine} 25: 4216-4226.

}
\author{Georg Heinze}
\seealso{
\code{\link{anova.logistf}}
}
\examples{
data(sex2)
fit<-logistf(case ~ age+oc+vic+vicl+vis+dia,  data=sex2)

logistftest(fit, test = ~ vic + vicl - 1, values = c(2, 0))
                        
}
\keyword{regression}
\keyword{models}
