rk <- 
  function(x, df = NULL, knots = NULL, m = NULL, intercept = FALSE, 
           Boundary.knots = NULL, warn.outside = TRUE, 
           periodic = FALSE, xlev = levels(x)){
    # reproducing kernel spline basis
    # Nathaniel E. Helwig (helwig@umn.edu)
    # 2023-09-06
    
    
    ######***######   NOMINAL BASIS   ######***######
    
    if(!is.null(xlev) && !is.ordered(x)){
      
      # check x
      x <- factor(x, levels = xlev)
      nlev <- length(xlev)
      n <- length(x)
      
      # check df and knots
      dropone <- FALSE
      if(is.null(knots)){
        if(is.null(df)) {
          df <- nlev
        } else {
          df <- as.integer(df[1])
          if(df < 1L) stop("The 'df' argument must satisfy:  df >= 1")
          if(df > nlev) stop("The 'df' argument must satisfy:  df <= nlevels(x)")
        }
        knots <- xlev[seq(1, nlev, length.out = df)]
        knots <- factor(knots, levels = xlev)
        dropone <- TRUE
      } else {
        knots <- factor(knots, levels = xlev)
        df <- length(knots)
        knotchar <- as.character(sort(knots))
        if((df == nlev) && identical(knotchar, xlev)) dropone <- TRUE
      }
      
      # evaluate basis
      mat <- outer(x, knots, FUN = "==") + 0.0
      colnames(mat) <- knots
      
      # drop last column?
      if(dropone){
        mat <- mat[,1:(nlev-1),drop=FALSE]
        idx <- which(x == knots[nlev])
        if(length(idx) > 0) mat[idx,] <- -1.0
      }
      
      # intercept?
      if(intercept) {
        cnames <- c("(Intercept)", paste0("knot", 1:ncol(mat)))
        mat <- cbind(1, mat)
      } else {
        cnames <- paste0("knot", 1:ncol(mat))
      }
      colnames(mat) <- cnames
      
      # append attributes
      attr(mat, "df") <- df
      attr(mat, "knots") <- knots
      attr(mat, "m") <- 0L
      attr(mat, "intercept") <- intercept
      attr(mat, "Boundary.knots") <- Boundary.knots
      attr(mat, "periodic") <- periodic
      attr(mat, "xlev") <- xlev
      
      # return basis
      return(mat)
      
    } # end if(!is.null(xlev) && !is.ordered(x))
    
    
    ######***######   ORDINAL BASIS   ######***######
    
    # check m
    if(is.null(m)){
      m <- ifelse(is.null(xlev), 2L, 0L)
    } else {
      m <- as.integer(m[1])
      if(m < 0L | m > 3L) stop("Input 'm' must be 0 (ordinal), 1 (linear), 2 (cubic), or 3 (quintic)")
    }
    
    # is ordinal?
    if(m == 0L){
      
      # check x
      if(is.null(xlev)){
        x <- as.ordered(x)
        xlev <- levels(x)
      } else {
        x <- factor(x, levels = xlev, ordered = TRUE)
      }
      nlev <- length(xlev)
      n <- length(x)
      
      # check df and knots
      if(is.null(knots)){
        if(is.null(df)) {
          df <- nlev
        } else {
          df <- as.integer(df[1])
          if(df < 1L) stop("The 'df' argument must satisfy:  df >= 1")
          if(df > nlev) stop("The 'df' argument must satisfy:  df <= nlevels(x)")
        }
        knots <- xlev[seq(1, nlev, length.out = df)]
        knots <- factor(knots, levels = xlev, ordered = TRUE)
      } else {
        knots <- factor(knots, levels = xlev, ordered = TRUE)
        df <- length(knots)
      }
      
      # transform x and knots
      ix <- as.integer(x)
      ik <- as.integer(knots)
      
      # evaluate basis
      const <- (nlev - 1) * (2*nlev - 1) / (6 * nlev)
      mat <- 1 - outer(ix, ik, pmax) + const
      mat <- mat + outer(ix*(ix-1), ik*(ik-1), FUN = "+") / (2 * nlev)
      
      # evaluate penalty
      pen <- 1 - outer(ik, ik, pmax) + const
      pen <- pen + outer(ik*(ik-1), ik*(ik-1), FUN = "+") / (2 * nlev)
      
      # absorb (inverse square root of) penalty into basis
      peig <- eigen(pen, symmetric = TRUE)
      prank <- sum(peig$values > peig$values[1] * df * .Machine$double.eps)
      pisqrt <- peig$vectors[,1:prank,drop=FALSE] %*% diag(1/sqrt(peig$values[1:prank]), nrow = prank, ncol = prank)
      mat <- mat %*% pisqrt
      
      # intercept?
      if(intercept) {
        cnames <- c("(Intercept)", paste0("knot", 1:ncol(mat)))
        mat <- cbind(1, mat)
      } else {
        cnames <- paste0("knot", 1:ncol(mat))
      }
      colnames(mat) <- cnames
      
      # append attributes
      attr(mat, "df") <- df
      attr(mat, "knots") <- knots
      attr(mat, "m") <- m
      attr(mat, "intercept") <- intercept
      attr(mat, "Boundary.knots") <- Boundary.knots
      attr(mat, "periodic") <- periodic
      attr(mat, "xlev") <- xlev
      
      # return basis
      return(mat)
      
    } # end if(m == 0L)
    
    
    
    ######***######   POLYNOMIAL BASIS   ######***######
    
    # check x
    x <- as.numeric(x)
    n <- length(x)
    
    # check Boundary.knots
    if(is.null(Boundary.knots)){
      Boundary.knots <- range(x)
    } else {
      Boundary.knots <- Boundary.knots[1:2]
      if(Boundary.knots[1] >= Boundary.knots[2]) stop("Input 'Boundary.knots' must satisfy: Boundary.knots[1] < Boundary.knots[2]")
      if(warn.outside){
        if(min(x) < Boundary.knots[1] - 1e-8) warning("min(x) < Boundary.knots[1]")
        if(max(x) > Boundary.knots[2] + 1e-8) warning("max(x) > Boundary.knots[2]")
      }
    }
    
    # check df and knots
    if(is.null(knots)){
      if(is.null(df)){
        df <- 5L
      } else {
        df <- as.integer(df[1])
        if(df < 1L) stop("The 'df' argument must satisfy:  df >= 1")
        if(df > n) warning("Requested 'df' is greater than the length of 'x'")
      }
      knots <- quantile(x, probs = seq(0, 1, length.out = df))
    } else {
      knots <- as.numeric(knots)
      knots <- unique(sort(c(Boundary.knots, knots)))
      df <- length(knots)
    }
    
    # check periodic
    periodic <- as.logical(periodic[1])
    
    # bernoulli polynomial function
    bernpoly <- function(x, m){
      val <- x - 1/2
      if(m == 2L){
        val <- (val^2 - 1/12) / 2
      } else if(m == 3L){
        val <- (val^3 - val / 4) / 6
      } else if(m == 4L){
        val <- (val^4 - val^2/2 + 7/240) / 24
      } else if(m == 5L){
        val <- (val^5 - (5/6) * val^3 - (7/48) * val) / 120
      } else if(m == 6L){
        val <- (val^6 - (5/4) * val^4 + (7/16) * val^2 - 31/1344) / 720
      }
      return(val)
    } # end bp
    
    # transform x and knots
    xr <- diff(Boundary.knots)
    x <- (x - Boundary.knots[1]) / xr
    knots <- (knots - Boundary.knots[1]) / xr
    
    # evaluate basis
    mat <- matrix(0, nrow = n, ncol = df)
    if(!periodic){
      for(k in 1:m){
        mat <- mat + outer(bernpoly(x, m = k), bernpoly(knots, m = k))
      }
    }
    dif <- abs(outer(x, knots, FUN = "-"))
    mat <- mat + (-1)^(m-1) * bernpoly(dif, m = 2*m)
    
    # evaluate penalty
    pen <- matrix(0, nrow = df, ncol = df)
    if(!periodic){
      for(k in 1:m){
        pen <- pen + outer(bernpoly(knots, m = k), bernpoly(knots, m = k))
      }
    }
    dif <- abs(outer(knots, knots, FUN = "-"))
    pen <- pen + (-1)^(m-1) * bernpoly(dif, m = 2*m)
    
    # absorb (inverse square root of) penalty into basis
    peig <- eigen(pen, symmetric = TRUE)
    prank <- sum(peig$values > peig$values[1] * df * .Machine$double.eps)
    pisqrt <- peig$vectors[,1:prank,drop=FALSE] %*% diag(1/sqrt(peig$values[1:prank]), nrow = prank, ncol = prank)
    mat <- mat %*% pisqrt
    
    # intercept?
    if(intercept) {
      cnames <- c("(Intercept)", paste0("knot", 1:ncol(mat)))
      mat <- cbind(1, mat)
    } else {
      cnames <- paste0("knot", 1:ncol(mat))
    }
    colnames(mat) <- cnames
    
    # append attributes
    attr(mat, "df") <- df
    attr(mat, "knots") <- knots * xr + Boundary.knots[1]
    attr(mat, "m") <- m
    attr(mat, "intercept") <- intercept
    attr(mat, "Boundary.knots") <- Boundary.knots
    attr(mat, "periodic") <- periodic
    
    # return basis
    return(mat)
    
  } # end rk