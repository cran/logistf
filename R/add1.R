#' Add or Drop All Possible Single Terms to/from a \code{logistf} Model
#'
#'Compute all the single terms in the scope argument that can be added to or dropped 
#'from the model, fit those models and compute a table of the changes in fit.
#'
#'\code{drop1} and \code{add1} generate a table where for each variable the penalized 
#'likelihood ratio chi-squared, the degrees of freedom, and the p-value for dropping/adding this variable are given.
#'
#' @param object A fitted \code{logistf, flic} or \code{flac} object
#' @param scope The scope of variables considered for adding or dropping. Should be a 
#' vector of variable names. Can be left missing; the method will then use all variables 
#' in the object's data slot which are not identified as the response variable.
#' @param data The data frame used to fit the object.
#' @param test The type of test statistic. Currently, only the PLR test (penalized likelihood 
#' ratio test) is allowed for logistf fits.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A matrix with \code{nvar} rows and 3 columns (Chisquared, degrees of freedom, p-value).
#' @export
#' 
#' @encoding UTF-8
#' @examples
#' data(sex2) 
#' fit<-logistf(data=sex2, case~1, pl=FALSE) 
#' add1(fit, scope=c("dia", "age"), data=sex2)
#'  
#' fit2<-logistf(data=sex2, case~age+oc+dia+vic+vicl+vis) 
#' drop1(fit2, data=sex2)
#' 
#' @importFrom stats add.scope
#' @importFrom stats update.formula
#' 
#' 
#' @rdname add1
#' @method add1 logistf
#' @exportS3Method add1 logistf
add1.logistf<-function(object, scope, data, test="PLR", ...){
  if(missing(scope)) stop("please provide scope: no terms in scope for adding to object")
  else if(is.numeric(scope)) scope<-attr(terms(object),"term.labels")[scope]
  else if(!is.character(scope)) scope <- add.scope(object, update.formula(object, scope))

  mf <- match.call(expand.dots =FALSE)
  m <- match(c("object", "scope", "test", "data"), names(mf), 0L)
  mf <- mf[c(1, m)]
  object <- eval(mf$object, parent.frame())
  
  whereisresponse<-(match(scope,colnames(model.frame(object))))
  if(sum(is.na(whereisresponse))==0) scope<-scope[-whereisresponse]
  #scope<-scope[is.na(match(scope, attr(terms(object),"term.labels")))]
  variables<-scope
  
  nvar<-length(variables)
  mat<-matrix(0,nvar,3)
  for(i in 1:nvar){
    newform<-as.formula(paste(object$formula,variables[i], sep="+"))
    fit2 <- update(object, formula=newform, data = data)
    res<-anova(object, fit2)
    mat[i,1]<-res$chisq
    mat[i,2]<-res$df
    mat[i,3]<-res$pval
  }
  rownames(mat)<-variables
  colnames(mat)<-c("ChiSq","df","P-value")
  return(mat)
}
#' @exportS3Method add1 flic
add1.flic<-function(object, scope, data, test="PLR", ...){
  add1.logistf(object, scope, data, ...)
}

#' @exportS3Method add1 flac
add1.flac<-function(object, scope, data, test="PLR", ...){
  add1.logistf(object, scope, data, ...)
}
#' @aliases drop1
#' @method drop1 logistf
#' @exportS3Method drop1 logistf
drop1.logistf<-function(object, scope, data, test="PLR", ...){
  mf <- match.call(expand.dots =FALSE)
  m <- match(c("object", "scope", "test"), names(mf), 0L)
  mf <- mf[c(1, m)]
  object <- eval(mf$object, parent.frame())
  
  extras <- list(...)
  if(!is.null(extras$full.penalty.vec)){
    full.penalty.vec <- extras$full.penalty.vec
  }
  else full.penalty.vec <- NULL
  if(missing(scope)){
    variables <-attr(terms(object),"term.labels")
  }
  else {
    variables<-scope
  }

  nvar<-length(variables)
  
  mat<-matrix(0,nvar,3) #initialise output
  for(i in 1:nvar){ #for every variable: omit from the object model fitin anova
    #full.penalty option of backward function
    if (!is.null(full.penalty.vec)&& nvar!=length(full.penalty.vec)){
      ind <- match(full.penalty.vec, attr(terms(object), "term.labels"))
      coltofit <- setdiff(0:length(attr(terms(object), "term.labels")), ind)+1 #determine which columns to fit for the full model
      newform <- as.formula(paste("~", paste(variables[i], paste(full.penalty.vec,collapse="+"), sep="+")))
      res<-anova(object, formula=newform, method="nested", col.fit.object=coltofit, ...)
    }
    else {
      newform<-as.formula(paste("~",variables[i]))
      res<-anova(object, formula=newform, method="nested",...)
    }
    mat[i,1]<-res$chisq
    mat[i,2]<-res$df
    mat[i,3]<-res$pval
  }
  rownames(mat)<-variables
  colnames(mat)<-c("ChiSq","df","P-value")
  return(mat)
}

#' @method drop1 flic
#' @exportS3Method drop1 flic
drop1.flic<-function(object, scope, data, test="PLR", ...){
  drop1.logistf(object, scope, data,  ...)
}

#' @method drop1 flac
#' @exportS3Method drop1 flac
drop1.flac<-function(object, scope, data, test="PLR", ...){
  drop1.logistf(object, scope, data, augmented_data=TRUE,...)
}