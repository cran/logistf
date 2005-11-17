.packageName <- "logistf"
#library(logistf)
#LOGISTF library by Meinhard Ploner, Daniela Dunkler, Harry Southworth, Georg Heinze, Medical University of Vienna
#any comments to georg.heinze@meduniwien.ac.at
#Version 1.03 (Build 2005.03.30)

logistf <-
function(formula = attr(data, "formula"), data = sys.parent(), pl = TRUE, alpha = 0.05,
    maxit = 25, maxhs = 5, epsilon = 0.0001, maxstep = 10, firth = TRUE, beta0)
{
    #n <- nrow(data)
    y <- as.vector(model.extract(model.frame(formula, data = data), response))
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix

    k <- ncol(x)    ## Anzahl Effekte

    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }


    beta <- c(log((sum(y)/n)/(1 - sum(y)/n)), rep(0, k - 1))
    if(!missing(beta0))
        beta[1] <- beta[1] - sum(x %*% beta0)/n
    iter <- 0
    pi <- as.vector(1/(1 + exp( - x %*% beta)))
    loglik <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))
    if(firth) {
        XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    }
    repeat {
        iter <- iter + 1
        XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        covs <- solve(Fisher)   ### (X' W  X) ^ -1
        H <- crossprod(XW2, covs) %*% XW2
        if(firth)
            U.star <- crossprod(x, y - pi + diag(H) * (0.5 - pi))
        else U.star <- crossprod(x, y - pi)
        delta <- as.vector(covs %*% U.star)
        mx <- max(abs(delta))/maxstep
        if(mx > 1)
            delta <- delta/mx
        beta <- beta + delta
        loglik.old <- loglik
        for(halfs in 1:maxhs) {
## Half-Steps
            pi <- as.vector(1/(1 + exp( - x %*% beta)))
            loglik <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))
            if(firth) {
                XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
                Fisher <- crossprod(t(XW2)) #### X' W  X
                loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
            }
            if(loglik > loglik.old)
                break
            beta <- beta - delta * 2^( - halfs)
    ##beta-Aenderung verkleinern
        }
        if(iter == maxit | sum(abs(delta)) <= epsilon)
            break
    }

    fit <- list(coefficients = beta, alpha = alpha, var = covs, df = (k-int), loglik =
        logistftest(formula, data, test=coltotest,firth=firth)$loglik, iter = iter, n = n, terms =
        colnames(x), y = y, formula = formula(formula), call=match.call())
    fit$linear.predictors <- as.vector(x %*% beta)
    fit$predict <- pi
    fit$hat.diag <- diag(H)
    if(firth)
        fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    vars <- diag(covs)
    if(pl) {
        LL.0 <- loglik - qchisq(1 - alpha, 1)/2
        fit$ci.lower <- fit$ci.upper <- rep(0, k)
        for(i in 1:k) {
            fit$ci.lower[i] <- logistpl(x, y, beta, i, LL.0, maxit, maxhs,
                epsilon, maxstep, firth, -1)$beta
            fit$ci.upper[i] <- logistpl(x, y, beta, i, LL.0, maxit, maxhs,
                epsilon, maxstep, firth, 1)$beta
            fit$prob[i] <- logistftest(formula, data, test = i, 0, maxit,
                maxhs, epsilon, maxstep, firth)$prob
        }
        fit$method.ci <- "Profile Likelihood"
    }
    else {
        fit$prob <- 1 - pchisq((beta^2/vars), 1)
        fit$method.ci <- "Wald"
        fit$ci.lower <- as.vector(beta + qnorm(alpha/2) * vars^0.5)
        fit$ci.upper <- as.vector(beta + qnorm(1 - alpha/2) * vars^0.5)
    }
    names(fit$prob) <- names(fit$ci.upper) <- names(fit$ci.lower) <- names(fit$
        coefficients) <- dimnames(covs)[[1]]
    attr(fit, "class") <- c("logistf")
    fit
}

####################################################################
################# logistfplot ######################################
#################             ######################################

logistfplot <- function(formula = attr(data, "formula"), data = sys.parent(), which, pitch = 0.05, limits,
                    alpha = 0.05, maxit = 25, maxhs = 5, epsilon = 0.0001, maxstep = 10, firth = TRUE,
                    legends = TRUE){

# by MP, 06.02.01
# which  ... righthand formula des zu plottenden Term (z.B. ~B oder ~A:D)
# pitch  ... distances between points in std's
# limits ... vector of MIN & MAX in std's, default=extremes of both CI's
#            +- 0.5 std. of beta
#

# Next line added by Harry Southworth, 22/10/02.
 if (missing(which)) stop("You must specify which (a one-sided formula).")

 fit <- logistf(formula = formula, data = data, alpha = alpha, maxit = maxit,
                maxhs = maxhs, epsilon = epsilon, maxstep = maxstep, firth = firth, pl = TRUE)
 coefs <- coef(fit) ## "normale" Koeffizienten
 covs <- fit$var ## Varianzen
# n <- nrow(data)
 resp <- model.extract(model.frame(formula, data = data), response)
 mm <- model.matrix(formula, data = data) ## Model-Matrix
 n <- nrow(mm)
 cov.name <- labels(mm)[[2]]
 k <- ncol(mm) #--> nun Berechnungen fuer Schleife
 cov.name2 <- labels(model.matrix(which, data = data))[[2]] ## Label des Test-Fakt.
 pos <- match(cov.name2, cov.name) ## Position des Testfakors
 std.pos <- diag(fit$var)[pos]^0.5
 if(missing(limits)) {
  lim.pl <- (c(fit$ci.lower[pos], fit$ci.upper[pos]) - coef(fit)[pos])/std.pos
  limits <- c(min(qnorm(alpha/2), lim.pl[1]) - 0.5, max(qnorm(1 - alpha/2), lim.pl[2]) + 0.5)
 }

 limits <- c(floor(limits[1]/pitch) * pitch, ceiling(limits[2]/pitch) * pitch)

 knots <- seq(limits[1], limits[2], pitch)
 nn <- length(knots)
 res <- matrix(knots, nn, 3) #initialisiere Werte
 dimnames(res) <- list(1:nn, c("std", cov.name2, "log-likelihood"))
 for(i in 1:nn) {
  res[i, 2] <- coefs[pos] + covs[pos, pos]^0.5 * knots[i]
  if(i == 1)
   xx <- logistftest(formula, data, test = which, values = res[i, 2], maxit=maxit,
            maxhs = maxhs, epsilon = epsilon, maxstep = maxstep, firth = firth)
  else xx <- logistftest(formula, data, test = which, values = res[i, 2], maxit = maxit,
            maxhs = maxhs, epsilon = epsilon, maxstep = maxstep, firth = firth, beta0 = xx$beta) ##verwende vorige Lsung!

  res[i, 3] <- xx$loglik[1]
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

 segments(coef(fit)[pos] - qnorm(alpha/2) * std.pos, yy[1], coef(fit)[pos] - qnorm(1 - alpha/2) *
            std.pos, yy[1], lty = 6) ##Wald-CI
 segments(fit$ci.lower[pos], yy[2], fit$ci.upper[pos], yy[2], lty = 8) ##prof.pen.lik.-CI

 axis(side = 3, at = res[ind, 2], labels = res[ind, 1])

 mtext(expression(paste("distance from ", hat(beta)," in multiples of ", hat(sigma))), side = 3, line = 3)
## mtext(expression(paste(beta, " of ", cov.name2)), side=1, line = 3)
 par(mai = act.par$mai)
 
 if (legends)
    {
     legend(x=coef(fit)[pos],
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
 if(fit$terms[1]!="(Intercept)")
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
function(formula = attr(data, "formula"), data = sys.parent(), test, values, maxit =
    25, maxhs = 5, epsilon = 0.0001, maxstep = 10, firth = TRUE, beta0)
{
    #n <- nrow(data)
    y <- model.extract(model.frame(formula, data = data), response)
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix
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

    beta <- c(log((sum(y)/n)/(1 - sum(y)/n)), rep(0, k - 1))
##berechne Startwerte
    iter <- 0
    loglik <- rep(0, 2)
    pi <- as.vector(1/(1 + exp( - x %*% beta)))
    if(missing(beta0)) {
################## coxphfplot braucht dies nicht! ###
#       loglik[2] <- sum(y * log(pi) + (1 - y) * log(1 - pi))
        loglik[2] <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))

        if(firth) {
            XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)
#### X' (W ^ 1/2)
            Fisher <- crossprod(t(XW2)) #### X' W  X
            loglik[2] <- loglik[2] + 0.5 * determinant(Fisher)$modulus[1]
        }
        repeat {
            iter <- iter + 1
            XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)
#### X' (W ^ 1/2)
            Fisher <- crossprod(t(XW2)) #### X' W  X
            covs <- solve(Fisher)   ### (X' W  X) ^ -1
            H <- crossprod(XW2, covs) %*% XW2
            if(firth)
                U.star <- crossprod(x, y - pi + diag(H) * (0.5 - pi))
            else U.star <- crossprod(x, y - pi)
            delta <- as.vector(covs %*% U.star)
            mx <- max(abs(delta))/maxstep
            if(mx > 1)
                delta <- delta/mx
            beta <- beta + delta
            loglik.old <- loglik[2]
            for(halfs in 1:maxhs) {
## 5 Half-Steps
                pi <- as.vector(1/(1 + exp( - x %*% beta)))
#                loglik[2] <- sum(y * log(pi) + (1 - y) * log(1 - pi))
                loglik[2] <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))

                if(firth) {
                  XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
                  Fisher <- crossprod(t(XW2))   #### X' W  X
                  loglik[2] <- loglik[2] + 0.5 * determinant(Fisher)$
                    modulus[1]
                }
                if(loglik[2] > loglik.old)
                  break
                beta <- beta - delta * 2^( - halfs)
    ##beta-Aenderung verkleinern
            }
            if(iter == maxit | sum(abs(delta)) <= epsilon)
                break
        }
    }
########################################################
## Labels der Test-Fakt.
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
        offset <- beta0
    else offset <- rep(0, k)    ## Vektor der fixierten Werte
    if(!missing(values))
        offset[pos] <- values
    beta <- offset  ########################################
    iter <- 0
    pi <- as.vector(1/(1 + exp( - x %*% beta)))
#    loglik[1] <- sum(y * log(pi) + (1 - y) * log(1 - pi))
    loglik[1] <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))

    if(firth) {
        XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        loglik[1] <- loglik[1] + 0.5 * determinant(Fisher)$modulus[1]
    }
    repeat {
        if(k2 == k) break   ## -> Overall Test
        iter <- iter + 1
        XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        covs <- solve(Fisher)   ### (X' W  X) ^ -1
        H <- crossprod(XW2, covs) %*% XW2
        if(firth)
            U.star <- crossprod(x, y - pi + diag(H) * (0.5 - pi))
        else U.star <- crossprod(x, y - pi)
        XX.XW2 <- crossprod(x[,  - pos, drop = FALSE], diag(pi * (1 - pi))^0.5)
    #### Teil von X' (W ^ 1/2)
        XX.Fisher <- crossprod(t(XX.XW2))   #### Teil von  X' W  X
        XX.covs <- matrix(0, k, k)
        XX.covs[ - pos,  - pos] <- solve(XX.Fisher)
    ### aufblasen der Cov-Matrix
        delta <- as.vector(XX.covs %*% U.star)
        mx <- max(abs(delta))/maxstep
        if(mx > 1)
            delta <- delta/mx
        beta <- beta + delta
        loglik.old <- loglik[1]
        for(halfs in 1:maxhs) {
## Half-Steps
            pi <- as.vector(1/(1 + exp( - x %*% beta)))
#            loglik[1] <- sum(y * log(pi) + (1 - y) * log(1 - pi))
            loglik[1] <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))

            if(firth) {
                XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
                Fisher <- crossprod(t(XW2)) #### X' W  X
                loglik[1] <- loglik[1] + 0.5 * determinant(Fisher)$
                  modulus[1]
            }
            if(loglik[1] > loglik.old)
                break
            beta <- beta - delta * 2^( - halfs)
    ##beta-Aenderung verkleinern
        }
        if(iter == maxit | sum(abs(delta)) <= epsilon)
            break
    }
#######################
    offset[ - pos] <- NA
    names(offset) <- cov.name
    fit <- list(testcov = offset, loglik = loglik, df = k2, prob = 1 - pchisq(2 *
        diff(loglik), k2), call = match.call(), beta = beta)
    if(firth)
        fit$method <- "Penalized ML"
    else fit$method <- "Standard ML"
    attr(fit, "class") <- "logistftest"
    fit
}


logistpl <- function(x, y, beta, i, LL.0, maxit, maxhs, epsilon,
    maxstep, firth, which = -1)
{
## which -1...left, +1...right
    k <- length(beta)
    iter <- 0
    pi <- as.vector(1/(1 + exp( - x %*% beta)))
    XW2 <- crossprod(x, diag(pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
    Fisher <- crossprod(t(XW2)) #### X' W  X
#    loglik <- sum(y * log(pi) + (1 - y) * log(1 - pi))
    loglik <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))
   if(firth)
        loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    repeat {
        iter <- iter + 1
        covs <- solve(Fisher)
        H <- crossprod(XW2, covs) %*% XW2
        if(firth)
            U.star <- crossprod(x, y - pi +
                diag(H) * (0.5 - pi))
        else U.star <- crossprod(x, y - pi)
        V.inv <-  - covs
        lambda <- which * ((2 * ((LL.0 - loglik
            ) + 0.5 * crossprod(U.star,
            V.inv) %*% U.star))/V.inv[i, i]
            )^0.5
        delta <-  - V.inv %*% (U.star + lambda *

            diag(k)[i,  ])
        mx <- max(abs(delta))/maxstep
        if(mx > 1)
            delta <- delta/mx
        beta <- beta + delta
        loglik.old <- loglik
        pi <- as.vector(1/(1 + exp( - x %*%
            beta)))
#        loglik <- sum(y * log(pi) + (1 - y) *  log(1 - pi))
        loglik <- sum(log(pi[y==1]))+sum(log(1-pi[y==0]))

        if(firth) {
            XW2 <- crossprod(x, diag(pi * (
                1 - pi))^0.5)
    #### X' (W ^ 1/2)
            Fisher <- crossprod(t(XW2))
    #### X' W  X
            loglik <- loglik + 0.5 *
                determinant(Fisher)$
                modulus[1]
        }
        if(iter == maxit | abs(loglik - LL.0) <=

            epsilon)
            break
    }
    list(beta = beta[i], LL = loglik)
}

.First.lib <- function(...)
    cat("LOGISTF library by Meinhard Ploner,\n                   Daniela Dunkler*,\n                   Harry Southworth,\n                   Georg Heinze*,\n *Medical University of Vienna\nVersion 1.05 (Build 2005.11.17)\nfor comments mailto:georg.heinze@meduniwien.ac.at\n")
