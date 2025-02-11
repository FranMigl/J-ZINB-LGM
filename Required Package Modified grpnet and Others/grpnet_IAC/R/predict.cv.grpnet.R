# Some code (for extracting predictions) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

predict.cv.grpnet <-
  function(object, 
           newx,
           newdata,
           s = c("lambda.min", "lambda.1se"),
           type = c("link", "response", "class", "terms", "importance", 
                    "coefficients", "nonzero", "groups", "ncoefs", "ngroups", "norm"),
           ...){
    # predict from a fit cv.grpnet object
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-09-06
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### check object
    if(!inherits(object, "cv.grpnet")) stop("Input 'object' must be of class 'cv.grpnet'")
    
    ### check newx / newdata
    sometypes <- c("link", "response", "class", "terms", "importance")
    itype <- sometypes[pmatch(as.character(type[1]), sometypes)]
    if(itype %in% c("link", "response", "class", "terms", "importance")){
      if(is.null(object$grpnet.fit$formula)){
        newx <- as.matrix(newx)
        ncoefs <- if(is.list(object$grpnet.fit$beta)) nrow(object$grpnet.fit$beta[[1]]) else nrow(object$grpnet.fit$beta)
        if(ncol(newx) != ncoefs) stop("Input 'newx' must satisfy:  ncol(newx) == nrow(object$beta)")
      } else {
        vnames <- rownames(attr(terms(object$grpnet.fit$formula), "factors"))[-1]
        check <- match(vnames, names(newdata))
        if(any(is.na(check))) stop("Input 'data' is missing variables included in model formula")
        if(!is.null(object$grpnet.fit$rk.args)){
          newx <- rk.model.matrix(object = object$grpnet.fit$formula, data = newdata,
                                  knots = object$grpnet.fit$rk.args$knots,
                                  Boundary.knots = lapply(object$grpnet.fit$rk.args$knots, function(x) x[c(1, length(x))]),
                                  m = object$grpnet.fit$rk.args$m, 
                                  periodic = object$grpnet.fit$rk.args$periodic,
                                  xlev = object$grpnet.fit$rk.args$xlev)
        } else {
          newx <- model.matrix(object = object$grpnet.fit$formula, data = newdata)
        }
        if(colnames(newx)[1] == "(Intercept)") newx <- newx[,-1]
        ncoefs <- if(is.list(object$grpnet.fit$beta)) nrow(object$grpnet.fit$beta[[1]]) else nrow(object$grpnet.fit$beta)
        if(ncol(newx) != ncoefs) stop("Input 'newdata' produced a design matrix of the wrong dimension\n (likely due to a factor level mismatch between 'data' and 'newdata')")
        object$grpnet.fit$formula <- NULL
      }
      nobs <- nrow(newx)
      newxnames <- rownames(newx)
      if(is.null(newxnames)) newxnames <- 1:nobs
    } else {
      newx <- NULL
    }
    
    ### check s
    if(is.character(s)){
      if(!any(s[1] == c("lambda.min", "lambda.1se"))) stop("Input 's' must be a character ('lambda.min' or 'lambda.1se') or a numeric vector")
      id <- object$index[ifelse(s[1] == "lambda.min", 1, 2)]
      s <- object$lambda[id]
    } else if(is.numeric(s)){
      if(any(s < 0)) stop("Inpust 's' must be a vector of non-negative numerics")
    } else {
      stop("Input 's' must be a character ('lambda.min' or 'lambda.1se') or a numeric vector")
    } # end if(is.character(s))
    
    ### return predictions
    predict(object$grpnet.fit, newx = newx, s = s, type = type)
    
  } # end predict.cv.grpnet