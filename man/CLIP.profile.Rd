\name{CLIP.profile}
\alias{CLIP.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine Profile Likelihoods from Imputed-Data Model Fits
}
\description{
This function uses CLIP (combination of likelihood profiles) to compute the pooled profile of the posterior after multiple imputation.
}
\usage{
CLIP.profile(obj = NULL, variable, data, which, firth = TRUE, weightvar, 
   control = logistf.control(), offset = NULL, from = NULL, to = NULL, 
   steps = 101, legacy = FALSE, keep = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{which}{
Alternatively to \code{variable}, the argument \code{which} allows to specify the variable to compute the profile for as righthand formula, e.g. \code{which=~X}.
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
   \item{offset}{
An optional offset variable.
}
  \item{from}{
Lowest value for the sequence of values for the regression coefficients for which the profile will be computed. Can be left blank. 
}
  \item{to}{
Highest value for the sequence of values for the regression coefficients for which the profile will be computed. Can be left blank. 
}
  \item{steps}{
Number of steps for the sequence of values for the regression coefficients for which the profile will be computed. 
 }
  \item{legacy}{
If TRUE, only R code will be used. Avoid.
}
  \item{keep}{
If TRUE, keeps the profiles for each imputed data sets in the output object.
}
}
\details{
While CLIP.confint iterates to find those values at which the CDF of the pooled posterior equals the confidence levels, CLIP.profile will evaluate the whole profile, which enables plotting and evaluating the skewness of
the combined and the completed-data profiles. The combined and completed-data profiles are available as cumulative distribution function (CDF) or in the scaling of relative profile likelihood (minus twice the likelihood ratio
statistic compared to the maximum). Using a \code{plot} method, the pooled posterior can also be displayed as a density.

}
\value{
An object of class \code{CLIP.profile} with items
\item{beta}{the values of the regression coefficient}
\item{cdf}{the cumulative distribution function of the posterior}
\item{profile}{the profile of the posterior}
\item{cdf.matrix}{An imputations x steps matrix with the values of the completed-data CDFs for each beta}
\item{profile.matrix}{An imputations x steps matrix with the values of the completed-data profiles for each beta}
\item{call}{the function call}
}
\references{
Heinze G, Ploner M, Beyea J (2013). Confidence intervals after multiple imputation: combining profile likelihood information from logistic regressions. Statistics in Medicine, to appear. 
}
\author{
Georg Heinze and Meinhard Ploner}
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

# CLIP profile
xprof<-CLIP.profile(obj=fit.list, variable="x", keep=TRUE)
plot(xprof)

#plot as CDF
plot(xprof, "cdf")

#plot as density 
plot(xprof, "density")
}

\keyword{ models }
\keyword{ regression }% __ONLY ONE__ keyword per line
