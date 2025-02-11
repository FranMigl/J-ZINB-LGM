# Some code (for extracting coefficients) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

coef.cv.grpnet <-
  function(object, 
           s = c("lambda.min", "lambda.1se"),
           ...){
    # coefficients from a fit cv.grpnet object
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-07-05
    
    
    ### check object
    if(!inherits(object, "cv.grpnet")) stop("Input 'object' must be of class 'cv.grpnet'")
    
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
    coef(object$grpnet.fit, s = s, ...)
    
  } # end coef.cv.grpnet