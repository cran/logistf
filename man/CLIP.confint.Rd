\name{CLIP.confint}
\alias{CLIP.confint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals after multiple imputation: Combination of Likelihood Profiles
}
\description{
This function implements the new combination of likelihood profiles (CLIP) method described in Heinze, Ploner and Beyea (2013). This method is useful for computing confidence intervals for parameters after multiple imputation of data sets,
if the normality assumption on parameter estimates and consequently the validity of applying Rubin`s rules (pooling of variances) is in doubt. It consists of combining the profile likelihoods into a posterior. The function \code{CLIP.confint}
searches for those values of a regression coefficient, at which the cumulative distribution function of the posterior is equal to the values specified in the argument \code{ci.level} (usually 0.025 and 0.975).
The search is performed using R`s \code{optimize} function.
}
\usage{
CLIP.confint(obj = NULL, variable = NULL, data, firth = TRUE, weightvar = NULL, 
   control = logistf.control(), ci.level = c(0.025, 0.975), pvalue = TRUE,
   offset = NULL,    bound.lo = NULL, bound.up = NULL, legacy = FALSE)
}
\arguments{
  \item{obj}{
Either a list of logistf fits (on multiple imputed data sets), or the result of analysis of a  \code{mice} (multiply imputed) object using \code{with.mids}.
}
  \item{variable}{
The variable of interest, for which confidence intervals should be computed. If missing, confidence intervals for all variables will be computed.
}
  \item{data}{
A list of data set corresponding to the model fits. Can be left blank if \code{obj} was obtained with the \code{dataout=TRUE} option or if \code{obj} was obtained by \code{mice}.
}
  \item{firth}{
If TRUE, applies the Firth correction. Should correspond to the entry in \code{obj}.
}
  \item{weightvar}{
An optional weighting variable for each observation.
}
  \item{control}{
control parameters for \code{logistf}, usually obtained by \code{logistf.control()}
}
  \item{ci.level}{
The two confidence levels for each tail of the posterior distribution.
}
  \item{pvalue}{
If TRUE, will also compute a P-value from the posterior.
}
  \item{offset}{
An optional offset variable
}
  \item{bound.lo}{
Bounds (vector of length 2) for the lower limit. Can be left blank. Use only if problems are encountered.
}
  \item{bound.up}{
Bounds (vector of length 2) for the upper limit. Can be left blank. Use only if problems are encountered.
}
  \item{legacy}{
If TRUE, will use pure R code for all model fitting. Can be slow. Not recommended.
}
}
\details{
For each confidence limit, this function performs a binary search to evaluate the combined posterior, which is obtained by first transforming the imputed-data likelihood profiles into cumulative distribution functions (CDFs),
and then averaging the CDFs to obtain the CDF of the posterior. Usually, the binary search manages to find the confidence intervals very quickly. The number of iterations (mean and maximum) will be supplied in the output object.
Further details on the method can be found in Heinze, Ploner and Beyea (2013).
}
\value{
An object of class \code{CLIP.confint}, with items

\item{variable}{the variable(s) which were analyzed}
\item{estimate}{the pooled estimate (average over imputations}
\item{ci}{the confidence interval(s)}
\item{pvalue}{the pvalue(s)}
\item{imputations}{the number of imputed data sets}
\item{ci.level}{the confidence level (input)}
\item{bound.lo}{The bounds used for finding the lower confidence limit; usually not of interest. May be useful for error-tracing.}
\item{bound.up}{The bounds used for finding the upper confidence limit.}
\item{iter}{The number of iterations (for each variable and each tail (lower or upper)).}
\item{call}{the call to \code{CLIP.confint}}
}
\references{
Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining profile likelihood information from logistic regressions. Statistics in Medicine, to appear. 
}
\author{
Georg Heinze and Meinhard Ploner
}
\examples{
#generate data set with NAs
freq=c(5,2,2,7,5,4)
y<-c(rep(1,freq[1]+freq[2]), rep(0,freq[3]+freq[4]), rep(1,freq[5]), rep(0,freq[6]))
x<-c(rep(1,freq[1]), rep(0,freq[2]), rep(1,freq[3]), rep(0,freq[4]), rep(NA,freq[5]),
   rep(NA,freq[6]))
toy<-data.frame(x=x,y=y)


# impute data set 5 times
set.seed(169)
toymi<-list(0)
for(i in 1:5){
  toymi[[i]]<-toy
  y1<-toymi[[i]]$y==1 & is.na(toymi[[i]]$x)
  y0<-toymi[[i]]$y==0 & is.na(toymi[[i]]$x)
  xnew1<-rbinom(sum(y1),1,freq[1]/(freq[1]+freq[2]))
  xnew0<-rbinom(sum(y0),1,freq[3]/(freq[3]+freq[4]))
  toymi[[i]]$x[y1==TRUE]<-xnew1
  toymi[[i]]$x[y0==TRUE]<-xnew0
}


# logistf analyses of each imputed data set
fit.list<-lapply(1:5, function(X) logistf(data=toymi[[X]], y~x, pl=TRUE, dataout=TRUE))

# CLIP confidence limits
CLIP.confint(obj=fit.list)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
