.packageName <- "logistf"
#library(logistf)
#LOGISTF library by Meinhard Ploner, Daniela Dunkler, Harry Southworth, Georg Heinze, Medical University of Vienna
#any comments to georg.heinze@meduniwien.ac.at
#Version 1.1 (Build 2010.04.15)



logistf <-
function(formula = attr(data, "formula"), data = sys.parent(), pl = TRUE, alpha = 0.05,
    control, plcontrol, firth = TRUE, init, weights, plconf=NULL, ...)
{
    #n <- nrow(data)
#    if (is.null(weights)) weights<-rep(1,nrow(data))
   call <- match.call()
   if(missing(control)) control<-logistf.control()
   if(pl==TRUE & missing(plcontrol)) plcontrol<-logistpl.control()

    mf <- match.call(expand.dots =FALSE)
    m <- match(c("formula", "data","weights", "na.action", 
        "offset"), names(mf), 0L)
 #   mf<-model.frame(formula, data=data, weights=weights)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    y <- model.response(mf)
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix 
    k <- ncol(x)    ## Anzahl Effekte
    cov.name <- labels(x)[[2]]
    weight <- as.vector(model.weights(mf)  )
    offset <- as.vector(model.offset(mf)   )
    if (is.null(offset)) offset<-rep(0,n)
    if (is.null(weight)) weight<-rep(1,n)

    if (missing(init)) init<-rep(0,k)
    if (is.null(plconf) & pl==TRUE) plconf<-1:k

    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }

    fit.full<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=1:k, init, control=control)
    fit.null<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=int, init, control=control)
    fit <- list(coefficients = fit.full$beta, alpha = alpha, var = fit.full$var, df = (k-int), loglik =c(fit.null$loglik, fit.full$loglik),
        iter = fit.full$iter, n = n, terms =
        colnames(x), y = y, formula = formula(formula), call=match.call(), conv=fit.full$conv)
    names(fit$conv)<-c("LL change","max abs score","beta change")
    beta<-fit.full$beta
    covs<-fit.full$var
    pi<-fit.full$pi
    fit$linear.predictors <- as.vector(x %*% beta + offset)
    fit$predict <- fit.full$pi
    fit$hat.diag <- fit.full$hat.diag
    if(firth)
        fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    vars <- diag(covs)
    if(pl) {
        betahist.lo<-vector(length(plconf),mode="list")
        betahist.up<-vector(length(plconf),mode="list")
        pl.conv<-matrix(0,length(plconf),4)
        dimnames(pl.conv)[[1]]<-as.list(plconf)
        dimnames(pl.conv)[[2]]<-as.list(c("lower, loglik","lower, beta", "upper, loglik", "upper, beta"))
        LL.0 <- fit.full$loglik - qchisq(1 - alpha, 1)/2
        pl.iter<-matrix(0,k,2)
        fit$ci.lower <- fit$ci.upper <- rep(0, k)
        icount<-0
        for(i in plconf) {
            icount<-icount+1
            inter<-logistpl(x, y, beta, i, LL.0, firth, -1, offset, weight, plcontrol)
            fit$ci.lower[i] <- inter$beta
            pl.iter[i,1]<-inter$iter
            betahist.lo[[icount]]<-inter$betahist
            pl.conv.lower<-t(inter$conv)
            inter<-logistpl(x, y, beta, i, LL.0, firth, 1, offset, weight, plcontrol)
            fit$ci.upper[i] <- inter$beta
            pl.iter[i,2]<-inter$iter
            betahist.up[[icount]]<-inter$betahist
            pl.conv.upper<-t(inter$conv)
            pl.conv[i,]<-cbind(pl.conv.lower,pl.conv.upper)
            fit.i<-logistf.fit(x,y, weight=weight, offset=offset, firth, col.fit=(1:k)[-i], control=control)
            fit$prob[i] <- 1-pchisq(2*(fit.full$loglik-fit.i$loglik),1)
        }
        fit$pl.iter<-pl.iter
        fit$method.ci <- "Profile Likelihood"
        fit$betahist<-list(lower=betahist.lo, upper=betahist.up)
        fit$pl.conv<-pl.conv
    }
    else {
        fit$prob <- 1 - pchisq((beta^2/vars), 1)
        fit$method.ci <- "Wald"
        fit$ci.lower <- as.vector(beta + qnorm(alpha/2) * vars^0.5)
        fit$ci.upper <- as.vector(beta + qnorm(1 - alpha/2) * vars^0.5)
    }
    names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- names(fit$
        coefficients) <- dimnames(x)[[2]]
    attr(fit, "class") <- c("logistf")
    fit
}

####################################################################
################# logistfplot ######################################
#################             ######################################

logistfplot <- function(formula = attr(data, "formula"), data = sys.parent(), which, pitch = 0.05, limits,
                    alpha = 0.05,  firth = TRUE,
                    legends = TRUE, weights, control, plcontrol){

# by MP, 06.02.01
# which  ... righthand formula des zu plottenden Term (z.B. ~B oder ~A:D)
# pitch  ... distances between points in std's
# limits ... vector of MIN & MAX in std's, default=extremes of both CI's
#            +- 0.5 std. of beta
#

# Next line added by Harry Southworth, 22/10/02.
 if (missing(which)) stop("You must specify which (a one-sided formula).")
 if (missing(control)) control<-logistf.control()
 if (missing(plcontrol)) plcontrol<-logistpl.control()
 
   call <- match.call()

    mf <- match.call(expand.dots =FALSE)
    m <- match(c("formula", "data","weights", "na.action", 
        "offset"), names(mf), 0L)
 #   mf<-model.frame(formula, data=data, weights=weights)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    y <- model.response(mf)
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix 
    cov.name <- labels(x)[[2]]
    weight <- as.vector(model.weights(mf)  )
    offset <- as.vector(model.offset(mf)   )
    if (is.null(offset)) offset<-rep(0,n)
    if (is.null(weight)) weight<-rep(1,n)

    
    cov.name <- labels(x)[[2]]
    k <- ncol(x)
    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }
  cov.name2 <- labels(model.matrix(which, data = data))[[2]] ## Label des Test-Fakt.
  pos <- match(cov.name2, cov.name) ## Position des Testfakors
  fit<-logistf.fit(x, y, weight=weight, offset=offset, firth=firth, control=control) 
  std.pos <- diag(fit$var)[pos]^0.5
 
 coefs <- fit$beta ## "normale" Koeffizienten
 covs <- fit$var ## Varianzen
# n <- nrow(data)
 n <- nrow(x)
 cov.name <- labels(x)[[2]]
 if(missing(limits)) {
  lim.pl<-numeric(0)
  LL.0 <- fit$loglik - qchisq(1 - alpha, 1)/2
  lower.fit<-logistpl(x, y, init=fit$beta, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=-1, i=pos, plcontrol=plcontrol)
  lim.pl[1]<-lower.fit$beta
  upper.fit<-logistpl(x, y, init=fit$beta, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=1, i=pos, plcontrol=plcontrol)
  lim.pl[2]<-upper.fit$beta
  lim.pl <- (lim.pl - coefs[pos])/std.pos
  limits <- c(min(qnorm(alpha/2), lim.pl[1]) - 0.5, max(qnorm(1 - alpha/2), lim.pl[2]) + 0.5)
 }

 limits <- c(floor(limits[1]/pitch) * pitch, ceiling(limits[2]/pitch) * pitch)

 knots <- seq(limits[1], limits[2], pitch)
 nn <- length(knots)
 res <- matrix(knots, nn, 3) #initialisiere Werte
 dimnames(res) <- list(1:nn, c("std", cov.name2, "log-likelihood"))
 for(i in 1:nn) {
  res[i, 2] <- coefs[pos] + covs[pos, pos]^0.5 * knots[i]
  if(i == 1){
     init<-lower.fit$betahist[nrow(lower.fit$betahist),]
     init[pos]<-res[i,2]
     xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, col.fit<-(1:k)[-pos], init=init,
                   control=control) 
  }     
  else {
     init<-xx$beta
     init[pos]<-res[i,2]
     xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, col.fit<-(1:k)[-pos], init=init,
                   control=control) # use solution from last step
  }
  res[i, 3] <- xx$loglik
 }

 #### Graphischer Output:

 my.par <- act.par <- par()
 my.par$mai[3] <- 1.65 * act.par$mai[3]
## if(legends) my.par$mai[1] <- 2 * act.par$mai[1]
 par(mai = my.par$mai)
 ind <- (1:nn)[round(4 * res[, 1]) == round(4 * res[, 1], 10)]
 if(length(ind) == 0) ind <- 1:nn
 pp <- max(res[, 3]) - 0.5 * res[, 1]^2

 plot(res[, -1], type = "l", xlab=expression(beta)) ##Profile likelihood

 #lines(res[,2], pp, lty=4)  #<<<Wald approximative profile lik. >>>

 points(res[res[, 1] == 0, 2], max(res[, 3])) ##Maximum of likelihood

 segments(min(res[, 2]), max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1),
                max(res[, 2]), max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1), lty = 3) ##refer.line

 yy <- par("usr")[4] - (par("usr")[4] - par("usr")[3]) * c(0.9, 0.95)

 segments(fit$beta[pos] - qnorm(alpha/2) * std.pos, yy[1], fit$beta[pos] - qnorm(1 - alpha/2) *
            std.pos, yy[1], lty = 6) ##Wald-CI
 segments(lower.fit$beta, yy[2], upper.fit$beta, yy[2], lty = 8) ##prof.pen.lik.-CI

 axis(side = 3, at = res[ind, 2], labels = res[ind, 1])

 mtext(expression(paste("distance from ", hat(beta)," in multiples of ", hat(sigma))), side = 3, line = 3)
## mtext(expression(paste(beta, " of ", cov.name2)), side=1, line = 3)
 par(mai = act.par$mai)
 
 if (legends)
    {
     legend(x=fit$beta[pos],
            y=min((min(res[,3])+max(res[,3]))/2,(max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1))),
        legend=c("Profile penalized likelihood",
                 paste(100 * (1 - alpha),"%-reference line"),
                 "Wald confidence interval",
                 "Profile likelihood confidence interval"),
        lty=c(1,3,6,8), 
        text.col=c("black","black","black","black"), ncol=1, bty="n", xjust=0.5)
    }


 title(paste("Profile of penalized likelihood for Variable",cov.name2))
 invisible(res)
}

print.logistf <-
function(x, ...)
{
# x ... object of class logistf
 print(x$call)
 cat("Model fitted by", x$method)
 cat("\nConfidence intervals and p-values by", x$method.ci, "\n\n")
 out <- cbind(x$coefficients, diag(x$var)^0.5, x$ci.lower,
x$ci.upper, qchisq(1 - x$
  prob, 1), x$prob)
 dimnames(out) <- list(names(x$coefficients), c("coef", "se(coef)",
paste(c("lower", "upper"),
  1 - x$alpha), "z", "p"))
 if(x$method.ci != "Wald")
  dimnames(out)[[2]][5] <- "Chisq"
 print(out)
 LL <- 2 * diff(x$loglik)
 cat("\nLikelihood ratio test=", LL, " on ", x$df, " df, p=", 1 -
pchisq(LL, x$df), ", n=",
  x$n, "\n\n", sep = "")
 invisible(x)
}

summary.logistf <- function(object, ...){
# object ... object of class logistf
 print(object$call)
 cat("\nModel fitted by", object$method)
 cat("\nConfidence intervals and p-values by", object$method.ci, "\n\n")
 out <- cbind(object$coefficients, diag(object$var)^0.5, object$ci.lower,
object$ci.upper, qchisq(1 - object$
  prob, 1), object$prob)
 dimnames(out) <- list(names(object$coefficients), c("coef", "se(coef)",
paste(c("lower", "upper"),
  1 - object$alpha), "z", "p"))
 if(object$method.ci != "Wald")
  dimnames(out)[[2]][5] <- "Chisq"
 print(out)
 LL <- 2 * diff(object$loglik)
 cat("\nLikelihood ratio test=", LL, " on ", object$df, " df, p=", 1 -
pchisq(LL, object$df), ", n=",
  object$n, sep = "")
 if(object$terms[1]!="(Intercept)")
  wald.z <- t(coef(object)) %*% solve(object$var) %*% coef(object)
 else
  wald.z <- t(coef(object)[2:(object$df+1)]) %*%
            solve(object$var[2:(object$df+1),2:(object$df+1)]) %*%
            coef(object)[2:(object$df+1)]
 cat("\nWald test =", wald.z, "on", object$df, "df, p =", 1 - pchisq(wald.z, object$df))
 cat("\n\nCovariance-Matrix:\n")
 print(object$var)
 invisible(object)
}

####################################################################
################# logistftest ######################################
#################             ######################################


print.logistftest <- function(x, ...){
# x ... object of class logistftest
 print(x$call)
 cat("Model fitted by", x$method, "\n\nFactors fixed as follows:\n")
 print(x$testcov)
 LL <- 2 * diff(x$loglik)
 out <- c(x$loglik[1], x$loglik[2], LL/2)
 names(out) <- c("Restricted model", "Full model", "difference")
 cat("\nLikelihoods:\n")
 print(out)
 cat("\nLikelihood ratio test=", LL, " on ", x$df, " df, p=", x$prob, "\n", sep = "")
 invisible(x)
}


logistftest <-
function(formula = attr(data, "formula"), data = sys.parent(), test, values, firth = TRUE, beta0, weights, control)
{
   call <- match.call()
    if (missing(control)) control<-logistf.control()
    mf <- match.call(expand.dots =FALSE)
    m <- match(c("formula", "data","weights", "na.action", 
        "offset"), names(mf), 0L)
 #   mf<-model.frame(formula, data=data, weights=weights)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    y <- model.response(mf)
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix 
    cov.name <- labels(x)[[2]]
    weight <- as.vector(model.weights(mf)  )
    offset <- as.vector(model.offset(mf)   )
    if (is.null(offset)) offset<-rep(0,n)
    if (is.null(weight)) weight<-rep(1,n)

    
    cov.name <- labels(x)[[2]]
    k <- ncol(x)
    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }

###    fit.full<-logistf.fit(    ) # unrestricted, define init and col.fit from values, beta0 and test
###    fit.null<-logistf.fit(    ) # restricted, define init and col.fit from values, beta0 and test
    
    fit.full<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=1:k, control=control)
    
    pos<-coltotest
    if(missing(test))
        test <- coltotest
    if(is.vector(test))
        cov.name2 <- cov.name[test]
    else cov.name2 <- labels(model.matrix(test, data = data))[[2]]
    pos <- match(cov.name2, cov.name)   ## Position der Testfakt.
    OK <- !is.na(pos)
    pos <- pos[OK]
    cov.name2 <- cov.name2[OK]
    k2 <- length(cov.name2) ## Anzahl Faktoren
    if(!missing(beta0))
        offset1 <- beta0
    else offset1 <- rep(0, k)    ## Vektor der fixierten Werte
    if(!missing(values))
        offset1[pos] <- values
    beta <- offset1  ########################################

    fit.null<-logistf.fit(x=x, y=y, weight=weight, offset=offset, firth, col.fit=(1:k)[-pos], control=control, init=beta)

    loglik<-c(fit.null$loglik,fit.full$loglik)
    
    offset1[ - pos] <- NA
    names(offset1) <- cov.name
    fit <- list(testcov = offset1, loglik = loglik, df = k2, prob = 1 - pchisq(2 *
        diff(loglik), k2), call = match.call(), beta = beta)
    if(firth)
        fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    attr(fit, "class") <- "logistftest"
    fit
}



.First.lib <- function(...)
    cat("LOGISTF library by Meinhard Ploner,\n                   Daniela Dunkler*,\n                   Harry Southworth,\n                   Georg Heinze*,\n *Medical University of Vienna\nVersion 1.10 (Build 2010.04.13)\nfor comments mailto:georg.heinze@meduniwien.ac.at\n")
