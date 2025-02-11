# Some code (for interpolating coefficients) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

coef.grpnet <-
  function(object, 
           s = NULL,
           ...){
    # predict from a fit grpnet object
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-07-05
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### check object
    if(!inherits(object, "grpnet")) stop("Input 'object' must be of class 'grpnet'")
    
    ### check s
    if(is.null(s)){
      newlambdas <- FALSE
      s <- object$lambda
    } else {
      newlambdas <- TRUE
      s <- as.numeric(s)
      if(any(s < 0)) stop("Inpust 's' must be a vector of non-negative numerics")
    } # end if(is.null(s))
    
    ### number of fit lambdas
    nlam <- length(object$lambda)
    
    ### number of new lambdas
    ns <- length(s)
    
    
    ######***######   MULTINOMIAL   ######***######
    
    if(object$family$family == "multinomial"){ 
      
      ### get number of response classes
      nresp <- length(object$ylev)
      
      ### coefs at object$lambda are easy...
      if(!newlambdas){
        coefs <- vector("list", nresp)
        names(coefs) <- object$ylev
        for(k in 1:nresp){
          if(nlam == 1L){
            coefs[[k]] <- matrix(c(object$a0[k], object$beta[[k]]), ncol = 1)
            rownames(coefs[[k]]) <- c("(Intercept)", names(object$beta[[k]]))
            colnames(coefs[[k]]) <- "s1"
          } else {
            coefs[[k]] <- rbind(object$a0[k,], object$beta[[k]])
            rownames(coefs[[k]]) <- c("(Intercept)", rownames(object$beta[[k]]))
            colnames(coefs[[k]]) <- paste0("s", 1:nlam)
          }
        }
        class(coefs) <- "coef.grpnet"
        return(coefs)
      }
      
      ### only 1 fit lambda?
      if(nlam == 1L){
        coefs <- vector("list", nresp)
        names(coefs) <- object$ylev
        for(k in 1:nresp){
          coefs[[k]] <- matrix(c(object$a0[k], object$beta[[k]]), nrow = length(object$beta[[k]]) + 1, ncol = ns)
          rownames(coefs[[k]]) <- c("(Intercept)", names(object$beta[[k]]))
          colnames(coefs[[k]]) <- rep("s1", ns)
        }
        class(coefs) <- "coef.grpnet"
        return(coefs)
      }
      
      ### min and max lambda from fit model
      lambda.max <- max(object$lambda)
      lambda.min <- min(object$lambda)
      
      ### transform lambda (and s) to the [0,1] interval
      sfrac <- (lambda.max - s) / (lambda.max - lambda.min)
      lambda <- (lambda.max - object$lambda) / (lambda.max - lambda.min)
      
      ### correct for out-of-range s values
      sfrac[sfrac > max(lambda)] <- max(lambda)
      sfrac[sfrac < min(lambda)] <- min(lambda)
      
      ### linear interpolation at sfrac
      interp <- approx(x = lambda, y = seq(lambda), xout = sfrac)$y
      
      ### L and R indices
      left <- floor(interp)
      right <- ceiling(interp)
      
      ### update sfrac w.r.t. L and R indices
      sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
      sfrac[left == right] <- 1
      sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
      
      ### interpolate...
      coefs <- vector("list", nresp)
      names(coefs) <- object$ylev
      for(k in 1:nresp){
        coefs[[k]] <- rbind(object$a0[k,], object$beta[[k]])
        coefs[[k]] <- coefs[[k]][, left, drop=FALSE] %*% diag(sfrac, nrow = ns, ncol = ns) + coefs[[k]][, right, drop=FALSE] %*% diag(1 - sfrac, nrow = ns, ncol = ns)
        rownames(coefs[[k]]) <- c("(Intercept)", rownames(object$beta[[k]]))
        colnames(coefs[[k]]) <- paste0("s", 1:ns)
      }
      class(coefs) <- "coef.grpnet"
      return(coefs)
      
    } # end if(family == "multinomial")
    
    
    
    ######***######   OTHER FAMILIES   ######***######
    
    ### coefficients at object$lambda are easy...
    if(!newlambdas){
      if(nlam == 1L){
        coefs <- matrix(c(object$a0, object$beta), ncol = 1)
        rownames(coefs) <- c("(Intercept)", names(object$beta))
        colnames(coefs) <- "s1"
      } else {
        coefs <- rbind(object$a0, object$beta)
        rownames(coefs) <- c("(Intercept)", rownames(object$beta))
        colnames(coefs) <- colnames(object$beta)
      }
      class(coefs) <- "coef.grpnet"
      return(coefs)
    }
    
    ### only 1 fit lambda?
    if(nlam == 1L){
      coefs <- matrix(c(object$a0, object$beta), nrow = length(object$beta) + 1, ncol = ns)
      rownames(coefs) <- c("(Intercept)", names(object$beta))
      colnames(coefs) <- rep("s1", ns)
      class(coefs) <- "coef.grpnet"
      return(coefs)
    }
    
    ### min and max lambda from fit model
    lambda.max <- max(object$lambda)
    lambda.min <- min(object$lambda)
    
    ### transform lambda (and s) to the [0,1] interval
    sfrac <- (lambda.max - s) / (lambda.max - lambda.min)
    lambda <- (lambda.max - object$lambda) / (lambda.max - lambda.min)
    
    ### correct for out-of-range s values
    sfrac[sfrac > max(lambda)] <- max(lambda)
    sfrac[sfrac < min(lambda)] <- min(lambda)
    
    ### linear interpolation at sfrac
    interp <- approx(x = lambda, y = seq(lambda), xout = sfrac)$y
    
    ### L and R indices
    left <- floor(interp)
    right <- ceiling(interp)
    
    ### update sfrac w.r.t. L and R indices
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
    sfrac[abs(lambda[left] - lambda[right]) < .Machine$double.eps] <- 1
    
    ### interpolate...
    coefs <- rbind(object$a0, object$beta)
    coefs <- coefs[, left, drop=FALSE] %*% diag(sfrac, nrow = ns, ncol = ns) + coefs[, right, drop=FALSE] %*% diag(1 - sfrac, nrow = ns, ncol = ns)
    rownames(coefs) <- c("(Intercept)", rownames(object$beta))
    colnames(coefs) <- paste0("s", 1:ns)
    class(coefs) <- "coef.grpnet"
    return(coefs)
    
  } # end coef.grpnet

print.coef.grpnet <-
  function(x, ...){
    ndigits <- getOption("digits")
    if(is.list(x)){
      for(j in 1:length(x)){
        x[[j]] <- formatC(x[[j]], 
                          format = "f",
                          digits = ndigits,
                          zero.print = paste(rep(c(".", " "), times = c(1,  ndigits)), collapse = ""),
                          preserve.width = "common")
      }
    } else {
      x <- formatC(x, 
                   format = "f",
                   digits = ndigits,
                   zero.print = paste(rep(c(".", " "), times = c(1,  ndigits)), collapse = ""),
                   preserve.width = "common")
    }
    print.default(x, right = TRUE, quote = FALSE)
  }