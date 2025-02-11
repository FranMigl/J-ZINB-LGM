# Some code (for arguments and outputs) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

grpnet.formula <-
  function(formula,
           data,
           use.rk = TRUE,
           family = c("gaussian", "binomial", "multinomial", "poisson", 
                      "negative.binomial", "Gamma", "inverse.gaussian"),
           weights = NULL,
           offset = NULL,
           alpha = 1,
           nlambda = 100,
           lambda.min.ratio = ifelse(nobs < nvars, 0.05, 1e-4),
           lambda = NULL,
           penalty.factor = NULL,
           penalty = c("LASSO", "MCP", "SCAD"),
           gamma = ifelse(penalty == "MCP", 3, 4),
           theta = 1,
           standardize = TRUE,
           thresh = 1e-04,
           maxit = 1e05,
           ...){
    # group elastic net regularized regression (formula)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-09-05
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### get call
    grpnet.call <- match.call()
    
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
    
    ### penalty.factor
    pf <- NULL
    if(!is.null(penalty.factor)){
      pf <- sqrt(gsize)
      names(pf) <- term.labels
      penalty.factor <- as.list(penalty.factor)
      penfacname <- names(penalty.factor)
      if(is.null(penfacname)) stop("Input 'penalty.factor' must be a named list.")
      for(k in 1:length(penalty.factor)){
        id <- match(penfacname[k], term.labels)
        if(is.na(id)) warning("Input 'penalty.factor' contains the term:   '", penfacname[k], "'\nwhich does not match any of the model terms")
        pf[id] <- penalty.factor[[k]]
      }
    } # if(!is.null(penalty.factor))
    
    
    ######***######   GRPNET.DEFAULT   ######***######
    
    ### fit regularization path
    res <- grpnet.default(x = x, 
                          y = y, 
                          group = group,
                          family = family,
                          weights = weights,
                          offset = offset,
                          alpha = alpha, 
                          nlambda = nlambda,
                          lambda.min.ratio = lambda.min.ratio,
                          lambda = lambda,
                          penalty.factor = pf,
                          penalty = penalty,
                          gamma = gamma,
                          theta = theta,
                          standardize = standardize,
                          intercept = intercept,
                          thresh = thresh,
                          maxit = maxit)
    
    
    ######***######   POST-PROCESSING   ######***######
    
    ### correct the call
    res$call <- grpnet.call
    
    ### add the (potential expanded) formula
    res$formula <- as.formula(paste0(yname, " ~ ", paste(term.labels, collapse = " + ")))
    
    ### add the term.labels
    if(intercept) term.labels <- c("(Intercept)", term.labels)
    res$term.labels <- term.labels
    
    ### add the rk.args
    res$rk.args <- rk.args
    
    ### return results
    return(res)
    
  } # end grpnet.formula