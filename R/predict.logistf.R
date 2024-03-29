#' Predict Method for logistf Fits
#'
#' Obtains predictions from a fitted \code{logistf} object. 
#' 
#' If \code{newdata} is omitted the predictions are based on the data used for the fit. 
#' 
#' 
#' @param object A fitted object of class \code{logistf}.
#'
#' @param newdata Optionally, a data frame in which to look for variables with which to predict. 
#' If omitted, the fitted linear predictors are used.  
#' @param type The type of prediction required. The default is on the scale of the linear predictors. 
#' The alternative \code{response} gives the predicted probabilities. Type \code{terms} returns a matrix with the fitted
#' values of each term in the formula on the linear predictor scale.
#' @param flic If \code{TRUE}(default = \code{FALSE}), predictions are computed with intercept correction.
#' @param se.fit  If \code{TRUE}(default = \code{FALSE}) standard errors are computed.
#' @param reference  A named vector of reference values for each variable for \code{type="terms"}.
#' @param na.action Function determining what should be done with missing values in newdata. The default is to predict NA.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A vector or matrix of predictions.
#'
#' @rdname predict.logistf
#' @exportS3Method predict logistf

predict.logistf <- function (object, newdata, type = c("link", "response", "terms"), flic=FALSE, se.fit = FALSE, reference, na.action = na.pass,...) 
{
  predict_terms <- function(object, model_matrix){
    m <- model.matrix(object$formula, object$model)
    if(is.null(model_matrix)){
      mm <- m
    } else {
      mm <- model_matrix
    }
      aa <- attr(mm, "assign")
      ll <- attr(terms(object), "term.labels")
      ll <- c("(Intercept)", ll)
      aaa <- factor(aa, labels = ll)
      asgn <- split(order(aa), aaa)
      asgn$"(Intercept)" <- NULL
      avx <- colMeans(m)
      beta <- object$coefficients
      termsconst <- sum(avx * beta)
      nterms <- length(asgn)
      predictor <- matrix(ncol = nterms, nrow = NROW(mm))
      dimnames(predictor) <- list(rownames(mm), names(asgn))
      X <- sweep(mm, 2L, avx, check.margin = FALSE)
      for (i in seq.int(1L, nterms, length.out = nterms)) {
          iipiv <- asgn[[i]]
          predictor[, i] <- X[, iipiv, drop = FALSE] %*% beta[iipiv]
      }
      if(se.fit){
        se <- matrix(ncol = nterms, nrow = NROW(mm))
        dimnames(se) <- list(rownames(mm), names(asgn))
        Terms <- delete.response(terms(object))
        m <- model.frame(Terms, reference)
        reference <- model.matrix(Terms, m)
        for(t in attr(terms(object), "term.labels")){
          ind <- asgn[[t]]
          diffs <- mm[,ind,drop=FALSE]-reference[,ind]
          se_t <- apply(diffs, 1, function(x){
            t(x) %*% object$var[ind,ind] %*% x
          })
          se_t <- sqrt(se_t)
          se[, t] <- se_t
        }
      }
      attr(predictor, "constant") <- termsconst
      if (se.fit) {
        return(list(predictor = predictor, se.fit = se))
      } else {
        return(predictor)
      }
  }
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL

  X <- model.matrix(object$formula, object$model)
  if(type == "terms" && missing(reference) && se.fit){
    orig <- eval(object$call$data)[1,]
    names_orig <- names(orig)
    ind <- sapply(as.list(names_orig), grepl, x=as.character(object$formula), fixed=TRUE)
    orig <- orig[, ind]
    factor_pos <- sapply(orig, is.factor)
    factor_names <- names(orig)[factor_pos]
    for(i in factor_names){
      orig[1, i] <- levels(orig[,i])[1]
    }
    reference <- numeric(sum(ind)-sum(factor_pos))
    orig[1,!factor_pos] <- reference
    reference <- orig
    
  } else if(!missing(reference) && se.fit){
    reference <- data.frame(t(reference))
  }
  if (missing(newdata)){#no data - return linear.predictors or response according to type
    if (flic) {
      #check if flic=TRUE was set in object
      if(!object$flic) {
        message("predict called with flic=TRUE but logistf-object was called with flic=FALSE: refitting model for predictions")
        object <- update(object, flic=TRUE)
      }
    }
    pred <- switch(type, link = object$linear.predictors, response = object$predict, terms = predict_terms(object))
    if(se.fit && type!="terms"){
      se <- apply(X, 1, function(x){
              t(x) %*% object$var %*% x
      })
      se <- sqrt(se)
      se <- unname(se)
      if(type == "response"){
        ci_lower <- pred - 1.96*se
        ci_upper <- pred + 1.96*se
      }
    }
  }
  else {
    Terms <- delete.response(terms(object))
    m <- model.frame(Terms, newdata, na.action = na.action)
    X <-  model.matrix(Terms, m)
    if (flic) {
      if(!object$flic) {
        message("predict called with flic=TRUE but logistf-object was called with flic=FALSE: refitting model for predictions")
        object <- update(object, flic=TRUE, pl=FALSE)
      }
    }
    pred <- switch(type, link = as.numeric(object$coefficients %*% t(X)) , 
                     response = as.numeric(1/(1+exp(-(object$coefficients %*% t(X))))), 
                     terms = predict_terms(object, X))
    if(se.fit && type!="terms"){
      se <- apply(X, 1, function(x){
              t(x) %*% object$var %*% x
      })
      se <- sqrt(se)
      se <- unname(se)
      if(type == "response"){
        ci_lower <- pred - 1.96*se
        ci_upper <- pred + 1.96*se
      }
    }
  }
  if (se.fit) {
    if (type=="terms"){
      names(pred$predictor) <- NULL
      names(pred$se.fit) <- NULL
      return(list(fit = pred$predictor, se.fit = pred$se.fit))
    } else if(type=="response"){
      return(list(fit = (pred), lower = (ci_lower), upper=(ci_upper)))
    } else {
      return(list(fit = (pred), se.fit = (se)))
    }
  } else {
    return((pred))
  }
}