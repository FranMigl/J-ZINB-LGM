row.kronecker <- 
  function(X, Y){
    ### row-wise kronecker product
    ### Nathaniel E. Helwig (helwig@umn.edu)
    ### Aug 6, 2023
    
    nx <- nrow(X)
    ny <- nrow(Y)
    if(nx != ny) stop("Inputs 'X' and 'Y' must have the same number of rows")
    px <- ncol(X)
    py <- ncol(Y)
    xnames <- colnames(X)
    if(is.null(xnames)) xnames <- paste0("X", 1:px)
    ynames <- colnames(Y)
    if(is.null(ynames)) ynames <- paste0("Y", 1:py)
    znames <- NULL
    Z <- matrix(0.0, nrow = nx, ncol = px * py)
    for(j in 1:px){
      index <- 1:py + (j - 1) * py
      Z[,index] <- X[,j] * Y
      znames <- c(znames, paste0(xnames[j], ":", ynames))
    }
    colnames(Z) <- znames
    return(Z)
  }