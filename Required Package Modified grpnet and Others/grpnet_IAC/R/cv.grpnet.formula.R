# Some code (for arguments and outputs) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

cv.grpnet.formula <-
  function(formula,
           data, 
           use.rk = TRUE,
           weights = NULL,
           offset = NULL,
           alpha = c(0.01, 0.25, 0.5, 0.75, 1),
           type.measure = NULL,
           nfolds = 10, 
           foldid = NULL, 
           same.lambda = FALSE,
           parallel = FALSE, 
           cluster = NULL, 
           verbose = interactive(), 
           ...){
    # k-fold cross-validation for grpnet (formula)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-09-05
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### get call
    cv.grpnet.call <- match.call()
    
    ### check formula
    formula <- as.formula(formula)
    charform <- as.character(formula)
    if(charform[1] != "~") stop("Input 'formula' must be of the form:  y ~ x")
    
    ### model frame (modified from R's lm() function)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    ## need stats:: for non-standard evaluation
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    ### model response
    yname <- colnames(mf)[1]
    y <- model.response(mf, "any")
    
    ### model matrix
    rk.args <- NULL
    if(use.rk){
      x <- rk.model.matrix(object = formula, data = mf, ...)
      rk.args <- list(knots = attr(x, "knots"),
                      m = attr(x, "m"),
                      periodic = attr(x, "periodic"),
                      xlev = attr(x, "xlev"))
    } else {
      x <- model.matrix(object = formula, data = mf)
    }
    
    ### model terms
    terms <- attr(mf, "terms")
    xvar.class <- attr(terms, "dataClasses")[-1]
    term.labels <- attr(terms, "term.labels")
    
    ### nobs and nvars
    nobs <- nrow(x)
    nvars <- ncol(x)
    nxvars <- length(xvar.class)
    
    ### group vector
    group <- attr(x, "assign")
    gsize <- table(group)
    ngrps <- length(gsize)
    
    ### intercept?
    intercept <- ifelse(group[1] == 0, TRUE, FALSE)
    if(intercept) {
      x <- x[,-1]
      group <- group[-1]
      gsize <- gsize[-1]
      ngrps <- ngrps - 1L
      nvars <- nvars - 1L
    }
    
    
    ######***######   CV.GRPNET.DEFAULT   ######***######
    
    ### collect extra arguments
    args <- list(...)
    
    ### penalty.factor given?
    if(is.null(args$penalty.factor)){
      
      res <- cv.grpnet.default(x = x, 
                               y = y,
                               group = group,
                               weights = weights,
                               offset = offset,
                               alpha = alpha,
                               type.measure = type.measure,
                               nfolds = nfolds, 
                               foldid = foldid,
                               same.lambda = same.lambda,
                               parallel = parallel, 
                               cluster = cluster, 
                               verbose = verbose,
                               intercept = intercept, 
                               ...)
      
    } else {
      
      # prepare args$penalty.factor for input to cv.grpnet.default
      pf <- sqrt(gsize)
      names(pf) <- term.labels
      args$penalty.factor <- as.list(args$penalty.factor)
      penfacname <- names(args$penalty.factor)
      if(is.null(penfacname)) stop("Input 'penalty.factor' must be a named list.")
      for(k in 1:length(args$penalty.factor)){
        id <- match(penfacname[k], term.labels)
        if(is.na(id)) warning("Input 'penalty.factor' contains the term:   '", penfacname[k], "'\nwhich does not match any of the model terms")
        pf[id] <- args$penalty.factor[[k]]
      }
      args$penalty.factor <- pf
      
      # combine all arguments
      iargs <- list(x = x, 
                    y = y,
                    group = group,
                    weights = weights,
                    offset = offset,
                    alpha = alpha,
                    type.measure = type.measure,
                    nfolds = nfolds, 
                    foldid = foldid,
                    same.lambda = same.lambda,
                    parallel = parallel, 
                    cluster = cluster, 
                    verbose = verbose,
                    intercept = intercept)
      allargs <- c(iargs, args)
      
      # call cv.grpnet.default
      res <- do.call(cv.grpnet.default, allargs)
      
    } # end if(is.null(args$penalty.factor))
    
    
    
    ######***######   POST-PROCESSING   ######***######
    
    ### correct the call
    res$call <- cv.grpnet.call
    
    ### add the (potential expanded) formula
    res$grpnet.fit$formula <- as.formula(paste0(yname, " ~ ", paste(term.labels, collapse = " + ")))
    
    ### add the term.labels
    if(intercept) term.labels <- c("(Intercept)", term.labels)
    res$grpnet.fit$term.labels <- term.labels
    
    ### add the rk.args
    res$grpnet.fit$rk.args <- rk.args
    
    ### return results
    return(res)
    
  } # end cv.grpnet.formula