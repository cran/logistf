#-------------------------------------------------------------------------------
logDet    <- function(x) 2*sum(log(diag(chol(x))));
#-------------------------------------------------------------------------------
invFisher <- function(x) chol2inv(chol(x));
#-------------------------------------------------------------------------------



logistf.fit<-function(x, y, weight=NULL, offset=NULL, firth=TRUE, col.fit=NULL, init=NULL, control){
# fitter function for logistf

# lconv = convergence criterion for log likelihood
# gconv = convergence criterion for score
# xconv = convergence criterion for parameter estimates
# pos = columns in x which will not be estimated (left at init value)

n<-nrow(x)
k<-ncol(x)

if (is.null(init)) init=rep(0,k)
if (is.null(col.fit)) col.fit=1:k
if (is.null(offset)) offset=rep(0,n)
if (is.null(weight)) weight=rep(1,n)
if (col.fit[1]==0) maxit<-0   #only evaluate likelihood and go back
if (missing(control)) control<-logistf.control()

maxit<-control$maxit
maxstep<-control$maxstep
maxhs<-control$maxhs
lconv<-control$lconv
gconv<-control$gconv
xconv<-control$xconv

#pos <- col.fit
beta <- init
l.change <- 5

    iter <- 0
    pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
#    loglik <- sum(y * log(pi) + (1 - y) * log(1 - pi))
    loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))

    if(firth) {
        XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    }
    evals<-1
    repeat {
#        if(k2 == k) break   ## -> Overall Test
        loglik.old <- loglik
        beta.old <- beta
        XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)    #### X' (W ^ 1/2)
        Fisher <- crossprod(t(XW2)) #### X' W  X
        covs <- invFisher(Fisher)   ### (X' W  X) ^ -1
        H <- crossprod(XW2, covs) %*% XW2
        if(firth) U.star <- crossprod(x, weight*(y - pi) + diag(H) * (0.5 - pi))
        else U.star <- crossprod(x, weight*(y - pi))
        XX.covs <- matrix(0, k, k)
        if (col.fit[1] != 0){
               XX.XW2 <- crossprod(x[, col.fit, drop = FALSE], diag(weight * pi * (1 - pi))^0.5)
             #### Teil von X' (W ^ 1/2)
               XX.Fisher <- crossprod(t(XX.XW2))   #### Teil von  X' W  X
               XX.covs[col.fit, col.fit] <- invFisher(XX.Fisher)
              ### aufblasen der Cov-Matrix
        }
        delta <- as.vector(XX.covs %*% U.star)
        delta[is.na(delta)]<-0
        mx <- max(abs(delta))/maxstep
        if(mx > 1)
            delta <- delta/mx
         evals<-evals+1


        if(maxit > 0){
        iter <- iter + 1
        beta <- beta + delta
        for(halfs in 1:maxhs) {
## Half-Steps
            pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
            loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))
            if(firth) {
                XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
                #### X' (W ^ 1/2)
                Fisher <- crossprod(t(XW2)) #### X' W  X
#                loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
                loglik <- loglik + 0.5 * logDet(Fisher)

            }
            evals<-evals+1
            l.change<-loglik-loglik.old

            if(loglik > loglik.old)
                break
            beta <- beta - delta * 2^( - halfs)
            ##beta-Aenderung verkleinern
         }
        }
        if(iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star)<gconv)) & (all(l.change<lconv))))
            break
    }

    var <- XX.covs



 list(beta=beta, var=var, pi=pi, hat.diag=diag(H), loglik=loglik, iter=iter, evals=evals, conv=c(l.change, max(abs(U.star)), max(abs(delta))))
}





logistpl <- function(x, y, init=NULL, i, LL.0, firth, which = -1, offset=rep(0, length(y)), weight=rep(1,length(y)), plcontrol)
{
## which -1...left, +1...right
    k <- ncol(x)
    if (is.null(init)) init<-rep(0,k)
        beta<-init
    if (missing(plcontrol)) plcontrol<-logistpl.control()
    ortho<-plcontrol$ortho
    pr<-plcontrol$pr
    maxit<-plcontrol$maxit
    maxstep<-plcontrol$maxstep
    maxhs<-plcontrol$maxhs
    xconv<-plcontrol$xconv
    lconv<-plcontrol$lconv
    
    
    if(ortho==TRUE & k>1 & i>1){
       thecol<-x[,i]
       others<-x[,-i]
       x[,i]<-lm(thecol~others-1,weights=weight)$residuals
    }
    if(pr==TRUE & (k-1 > 1)){
       others<-x[,c(-1,-i)]
       pc1<-prcomp(others)
       x[,c(-1,-i)]<-predict(pc1,others)
    }
    
    iter <- 0
    betahist<-numeric(0)
    pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
    XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
    Fisher <- crossprod(t(XW2)) #### X' W  X
#    loglik <- sum(y * log(pi) + (1 - y) * log(1 - pi))
    loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))
   if(firth)
        loglik <- loglik + 0.5 * logDet(Fisher)
    repeat {
        iter <- iter + 1
        covs <- invFisher(Fisher)
        H <- crossprod(XW2, covs) %*% XW2
        if(firth)
            U.star <- crossprod(x, weight*(y - pi) +
                diag(H) * (0.5 - pi))
        else U.star <- crossprod(x, weight*(y - pi))
        V.inv <-  - covs
        lambda <- which * ((2 * ((LL.0 - loglik
            ) + 0.5 * crossprod(U.star,
            V.inv) %*% U.star))/V.inv[i, i]
            )^0.5
        delta <-  - V.inv %*% (U.star + lambda * diag(k)[i,  ])
        delta[is.na(delta)]<-0
        mx <- max(abs(delta))/maxstep
        if(mx > 1)
            delta <- delta/mx
        beta <- beta + delta
        loglik.old <- loglik
        hs<-0
        repeat {
         pi <- as.vector(1/(1 + exp( - x %*% beta - offset)))
#        loglik <- sum(y * log(pi) + (1 - y) *  log(1 - pi))
         loglik <- sum(weight[y==1]*log(pi[y==1]))+sum(weight[y==0]*log(1-pi[y==0]))

         if(firth) {
             XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
    #### X' (W ^ 1/2)
             Fisher <- crossprod(t(XW2))
    #### X' W  X
             loglik <- loglik + 0.5 * logDet(Fisher)
          }
         if((hs>maxhs)|((abs(loglik-LL.0)<abs(loglik.old-LL.0)) & (loglik>LL.0))) break
         beta<-beta - delta/2   ### a simple emergency step halfing
         delta<-delta/2
         hs<-hs+1
         }
        betahist<-rbind(betahist,t(beta))
        if(iter == maxit | ((abs(loglik - LL.0) <=  lconv)  & (max(abs(delta))<xconv)))
            break
    }
    list(beta = beta[i], LL = loglik, conv=c(abs(loglik - LL.0),  max(abs(delta))), iter=iter, betahist=betahist)
}

logistf.control<-function(maxit=25, maxhs=5, maxstep=5, lconv=0.00001, gconv=0.00001, xconv=0.00001){
  list(maxit=maxit, maxhs=maxhs, maxstep=maxstep, lconv=lconv, gconv=gconv, xconv=xconv)
}

logistpl.control<-function(maxit=100, maxhs=5, maxstep=5, lconv=0.00001, xconv=0.00001, ortho=FALSE, pr=FALSE){
  list(maxit=maxit, maxhs=maxhs, maxstep=maxstep, lconv=lconv, xconv=xconv, ortho=ortho, pr=pr)
}

