\name{pool.RR}
\Rdversion{1.1}
\alias{pool.RR}
\title{Compute Pooled Normal Confidence Intervals (following Rubin`s Rules) after Multiple Imputation}
\description{Computes pooled normal (symmetric) confidence intervals after multiple imputation, following Rubin`s rules.}
\usage{pool.RR(object, method = "plain")}
\arguments{
  \item{object}{A list of fits from logistf.}
  \item{method}{A string describing the method to compute the degrees of freedom, 
                \code{"plain"} produces the conventional degrees of freedom as in Rubin (1987), and this is required for fits from logistf.}
}
\details{
 This is a slightly modified version of the \code{pool()} function of the \code{mice} package, for use with a list of \code{logistf}-fits. 
 If logistf fits on multiple imputed versions of a data set are on hand in a list object, use this function to pool results using Rubin's rules.
 If logistf fits on multiple imputed versions of a data set are on hand in a \code{mira} object, use \code{mice}'s  \code{pool(..., method="plain")} function.

}
\value{
An object of class mipo, which stands for 'multiple imputation pooled outcome'. The object is a list containing the following items:
\item{call}{ 	The call to the pool function.}
\item{call1}{ 	The original call how the mira object was calculated.}
\item{call2}{ 	The original call to the mice function to calculate the underlying midsobject. (NA)}
\item{formula}{ 	The formula that is used in call1.}
\item{nmis}{ 	The number of missing entries for each variable. (NA)}
\item{m}{ 	The number of imputations}
\item{qhat}{ 	A matrix, containing the estimated coeffients of the m repeated complete data analyses}
\item{u}{ 	The corresponding m variancematrices of the estimates in an three dimensional array.}
\item{qbar}{ 	The pooled estimate, formula (3.1.2) Rubin (1987).}
\item{ubar}{ 	The mean of the variances, formula (3.1.3), Rubin (1987).}
\item{b}{ 	The within imputation variance, formula (3.1.4), Rubin (1987).}
\item{t}{ 	Total variance of the pooled estimates, formula (3.1.5), Rubin (1987).}
\item{r}{ 	Relative increase in variance due to nonresponse, formula (3.1.7), Rubin (1987).}
\item{df}{ 	Degrees of freedom for t reference distribution, calculated according to the article of Barnard and Rubin (1999).}
\item{f}{ 	Fraction missing information due to nonresponse, formula (3.1.10), Rubin (1987).}
}
\references{
Barnard, J. and Rubin, D.B. (1999). Small sample degrees of freedom with multiple imputation. Biometrika, 86, 948-955.

Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys. New York: John Wiley and Sons.

Van Buuren, S., Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation by Chained Equations in R. 
 Journal of Statistical Software, Volume 45, Issue 3.

Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS. Berlin: Springer. 

}
\author{adapted from \code{\link{mice}} package for use with \code{logistf} objects by Georg Heinze}

\seealso{
\code{\link{mice}}, \code{\link{CLIP.confint}}, \code{\link{PVR.confint}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
summary(pool.RR(fit.list))

}
\keyword{regression}
\keyword{models}
