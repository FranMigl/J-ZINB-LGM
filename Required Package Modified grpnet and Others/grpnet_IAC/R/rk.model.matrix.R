rk.model.matrix <- 
  function(object, data = environment(object), ...){
    # reproducing kernel model matrix
    # Nathaniel E. Helwig (helwig@umn.edu)
    # 2023-09-06
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### check object
    oc <- class(object)[1]
    if(!any(oc == c("formula", "terms"))){
      stop("Input 'object' should be a model formula or a terms object")
    }
    
    ### make terms
    if(oc == "formula"){
      formula <- object
      object <- terms(object, data = data)
    } else {
      formula <- object
      attributes(formula) <- NULL
      formula <- as.formula(formula)
    }
    
    ### extract information
    factors <- attr(object, "factors")
    nvars <- nrow(factors)
    nterms <- ncol(factors)
    variables <- rownames(factors)
    term.labels <- attr(object, "term.labels")
    order <- attr(object, "order")
    intercept <- attr(object, "intercept")
    response <- attr(object, "response")
    
    ### was response included?
    if(response == 1){
      variables <- variables[-1]
      nvars <- nvars - 1L
      factors <- factors[-1,,drop=FALSE]
    } else if(response > 0){
      stop("Input 'object' formula should be of the form:\n y ~ . (one response)\n   ~ . (no response)")
    }
    
    
    ######***######   ADDITIONAL ARGUMENTS   ######***######
    
    ### collect arguments
    args <- list(...)
    
    ### was df provided?
    if(is.null(args$df)){
      args$df <- vector("list", nvars)
      names(args$df) <- variables
    } else {
      if(length(args$df) == 1L){
        df <- args$df
        args$df <- vector("list", nvars)
        names(args$df) <- variables
        args$df[1:nvars] <- df
        rm(df)
      } else {
        if(!is.list(args$df)) stop("Input 'df' should be a single integer (common df for all variables)\n or a named list giving the df for each variable.")
      }
    } # end if(is.null(args$df))
    
    ### were knots provided?
    if(is.null(args$knots)){
      args$knots <- vector("list", nvars)
      names(args$knots) <- variables
    } else {
      if(!is.list(args$knots) | is.null(names(args$knots))) {
        stop("Input 'knots' should be a named list giving the knots for each variable.")
      }
    }
    
    ### was m provided?
    if(is.null(args$m)){
      args$m <- vector("list", nvars)
      names(args$m) <- variables
    } else {
      if(length(args$m) == 1L){
        m <- args$m
        args$m <- vector("list", nvars)
        names(args$m) <- variables
        args$m[1:nvars] <- m
        rm(m)
      } else {
        if(!is.list(args$m)) stop("Input 'm' should be a single integer (common m for all variables)\n or a named list giving the penalty order for each variable.")
      }
    } # end if(is.null(args$m))
    
    ### were Boundary.knots provided?
    if(is.null(args$Boundary.knots)){
      args$Boundary.knots <- vector("list", nvars)
      names(args$Boundary.knots) <- variables
    } else {
      if(!is.list(args$Boundary.knots) | is.null(names(args$Boundary.knots))) {
        stop("Input 'Boundary.knots' should be a named list giving the boundary knots for each variable.")
      }
    }
    
    ### was periodic provided?
    if(is.null(args$periodic)){
      args$periodic <- vector("list", nvars)
      args$periodic[1:nvars] <- FALSE
      names(args$periodic) <- variables
    } else {
      if(length(args$periodic) == 1L){
        periodic <- args$periodic
        args$periodic <- vector("list", nvars)
        names(args$periodic) <- variables
        args$periodic[1:nvars] <- periodic
        rm(periodic)
      } else {
        if(!is.list(args$periodic)) stop("Input 'periodic' should be a single logical (common periodicity info for all variables)\n or a named list giving the periodicity info for each variable.")
      }
    } # end if(is.null(args$periodic))
    
    ### was xlev provided?
    if(is.null(args$xlev)){
      args$xlev <- vector("list", nvars)
      names(args$xlev) <- variables
      for(j in 1:nvars){
        temp <- levels(data[,variables[j]])
        if(!is.null(temp)) args$xlev[[j]] <- temp
      }
    } else {
      if(!is.list(args$xlev) | is.null(names(args$xlev))) {
        stop("Input 'xlev' should be a named list giving the levels for each variable.")
      }
    }
    
    
    ######***######   BUILD MODEL MATRIX   ######***######
    
    ### make marginal basis matrices
    xlist <- vector("list", nvars)
    names(xlist) <- variables
    for(j in 1:nvars){
      xlist[[j]] <- rk(x = data[,variables[j]],
                       df = args$df[[variables[j]]],
                       knots = args$knots[[variables[j]]],
                       m = args$m[[variables[j]]],
                       Boundary.knots = args$Boundary.knots[[variables[j]]],
                       periodic = args$periodic[[variables[j]]],
                       xlev = args$xlev[[variables[j]]])
      colnames(xlist[[j]]) <- paste(variables[j], colnames(xlist[[j]]), sep = ".")
    }
    df.marg <- sapply(xlist, ncol)
    nobs <- nrow(xlist[[1]])
    
    ### find number of columns of x (and g)
    df.total <- 0
    g <- NULL
    for(k in 1:nterms){
      id <- which(factors[,k] > 0)
      df.new <- prod(df.marg[id])
      df.total <- df.total + df.new
      g <- c(g, rep(k, df.new))
    }
    
    ### build model matrix
    x <- matrix(0, nrow = nobs, ncol = df.total)
    xnames <- NULL
    for(k in 1:nterms){
      id <- which(g == k)
      if(order[k] == 1L){
        x[,id] <- xlist[[term.labels[k]]]
        xnames <- c(xnames, colnames(xlist[[term.labels[k]]]))
      } else {
        tid <- which(factors[,k] > 0)
        xtemp <- xlist[[variables[tid[1]]]]
        for(j in 2:order[k]){
          xtemp <- row.kronecker(xtemp, xlist[[variables[tid[j]]]])
        }
        xnames <- c(xnames, colnames(xtemp))
        x[,id] <- xtemp
      } # end if(order[k] == 1L)
    } # end for(k in 1:nterms)
    colnames(x) <- xnames
    
    ### add intercept?
    if(intercept > 0){
      x <- cbind(1, x)
      colnames(x) <- c("(Intercept)", xnames)
      g <- c(0, g)
    }
    
    
    ######***######   RETURN RESULTS   ######***######
    
    ### attributes
    attr(x, "assign") <- g
    attr(x, "term.labels") <- term.labels
    attr(x, "knots") <- lapply(xlist, function(x) attr(x, "knots"))
    attr(x, "m") <- lapply(xlist, function(x) attr(x, "m"))
    attr(x, "periodic") <- lapply(xlist, function(x) attr(x, "periodic"))
    attr(x, "xlev") <- lapply(xlist, function(x) attr(x, "xlev"))
    
    ### return
    return(x)
    
  } # end rk.model.matrix