\name{anova.logistf}
\alias{anova.logistf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analysis of Penalized Deviance for \code{logistf} Models
}
\description{
This method compares hierarchical and non-hierarchical logistf models using penalized likelhood ratio tests. It replaces the function logistftest of former versions of logistf.
}
\usage{
\method{anova}{logistf}(object, fit2, formula, method = "nested", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted logistf model object
}
  \item{fit2}{
another fitted logistf model object, to be compared with \code{object}
}
  \item{formula}{
alternatively to \code{fit2}, a formula which specifies terms to omit from the \code{object} model fit.
}
  \item{method}{
One of \code{c("nested","PLR")}. \code{nested} is the default for hierarchically nested models, and will compare the penalized likelihood ratio statistics (minus twice the difference between maximized penalized
log likelihood and null penalized log likelihood), where the null penalized log likelihood is computed from the same, hierarchically superior model. Note that unlike in maximum likelihood analysis, 
the null penalized likelihood depends on the penalty (Jeffreys prior) which itself depends on the scope of variables of the hierarchically superior model. \code{PLR} compares the difference in penalized likelihood ratio
between the two models, where for each model the null penalized likelihood is computed within the scope of variables in that model. For PLR, the models need not be hierarchically nested.
}
  \item{\dots}{
Further arguments passed to the method.
}
}
\details{
Comparing models fitted by penalized methods, one must consider that the penalized likelihoods are not directly comparable, since a penalty is involved. Or in other words, inserting zero for some regression coefficients
will not lead to the same penalized likelihood as if the corresponding variables are simply "unknown" to a model. The \code{anova} method takes care that the same penalty is used 
for two hierarchically nested models, and if the models are not hierarchically nested, it will first relate each penalized likelihood to its null penalized likelihood, and only compare the resulting penalized likelihod ratio
statistics. The chi-squared approximation for this latter method (\code{PLR}) is considered less accurate than that of the \code{nested} method. Nevertheless, it is the only way to go for comparison of non-nested models.


}
\value{
An object of class \code{anova.logistf} with items
\item{chisq}{the chisquared statistic for the model comparison}
\item{df}{the degrees of freedom}
\item{pval}{the p-value}
\item{call}{the function call}
\item{method}{the method of comparison (input)}
\item{model1}{the first model}
\item{model2}{the second model which was compared to the first model}
\item{PLR1}{the PLR statistic of the first model}
\item{PLR2}{the PLR statistic of the second model; for the \code{nested} method, this will be the drop in chi-squared due to setting the coefficients to zero}

}
\author{
Georg Heinze
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(sex2)
fit<-logistf(data=sex2, case~age+oc+dia+vic+vicl+vis)

#simultaneous test of variables vic, vicl, vis:
anova(fit, formula=~vic+vicl+vis)

#test versus a simpler model
fit2<-logistf(data=sex2, case~age+oc+dia)
# or: fit2<-update(fit, case~age+oc+dia)
anova(fit,fit2)


# comparison of non-nested models (with different df):
fit3<-logistf(data=sex2, case~age+vic+vicl+vis)
anova(fit2,fit3, method="PLR")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
