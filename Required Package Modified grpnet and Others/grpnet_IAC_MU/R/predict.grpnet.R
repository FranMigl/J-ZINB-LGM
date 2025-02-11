# Some code (for interpolating predictions) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

predict.grpnet <-
  function(object, 
           newx,
           newdata,
           s = NULL,
           type = c("link", "response", "class", "terms", "importance", 
                    "coefficients", "nonzero", "groups", "ncoefs", "ngroups", "norm"),
           ...){
    # predict from a fit grpnet object
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-09-06
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### check object
    if(!inherits(object, "grpnet")) stop("Input 'object' must be of class 'grpnet'")
  
    ### check type
    thetypes <- c("link", "response", "class", "terms", "importance", "coefficients", "nonzero", "groups", "ncoefs", "ngroups", "norm")
    family <- object$family$family
    type <- pmatch(as.character(type[1]), thetypes)
    if(is.na(type)) stop("Invalid 'type' input")
    type <- thetypes[type]
    if(type == "class" && !(family %in% c("binomial", "multinomial")))
      stop("Input 'type' can only be set to 'class' for binomial and multinomial families")
    
    
    
    ######***######   COEFFICIENTS & RELATED INFO   ######***######
    
    ### type == "coefficients"
    if(type == "coefficients"){
      return(coef.grpnet(object = object, s = s))
    }
    
    ### type == "nonzero"
    if(type == "nonzero"){
      coefs <- coef.grpnet(object = object, s = s)
      if(family == "multinomial"){
        return(apply(Reduce("+", lapply(coefs, function(x) abs(x) > 0)), 2, function(x) which(x > 0)))
      } else {
        return(apply(abs(coefs) > 0, 2, which))
      }
    }
    
    ### type == "ncoefs"
    if(type == "ncoefs"){
      coefs <- coef.grpnet(object = object, s = s)
      if(family == "multinomial"){
        return(Reduce("+", lapply(coefs, function(x) colSums(abs(x) > 0))))
      } else {
        return(colSums(abs(coefs) > 0))
      }
    }
    
    ### l2norm
    if(type %in% c("groups", "ngroups", "norm")){
      if(family == "multinomial"){
        l2norm <- function(x) sum(x^2)
      } else {
        l2norm <- function(x) sqrt(sum(x^2))
      }
      grpnorm <- function(x) tapply(x, object$group, l2norm)
    }
    
    ### type == "groups"
    if(type == "groups"){
      coefs <- coef.grpnet(object = object, s = s)
      if(family == "multinomial"){
        norms <- sqrt(Reduce("+", lapply(coefs, function(x) apply(x, 2, grpnorm))))
      } else {
        norms <- apply(coefs, 2, grpnorm)
      }
      if(ncol(norms) > 1L){
        return(apply(norms, 2, function(x) object$term.labels[which(x > 0)]))
      } else {
        return(object$term.labels[which(norms > 0)])
      }
    }
    
    ### type == "ngroups"
    if(type == "ngroups"){
      coefs <- coef.grpnet(object = object, s = s)
      if(family == "multinomial"){
        return(colSums(sqrt(Reduce("+", lapply(coefs, function(x) apply(x, 2, grpnorm)))) > 0))
      } else {
        return(colSums(apply(coefs, 2, grpnorm) > 0))
      }
    }
    
    ### type == "norm"
    if(type == "norm"){
      coefs <- coef.grpnet(object = object, s = s)
      if(family == "multinomial"){
        norms <- sqrt(Reduce("+", lapply(coefs, function(x) apply(x, 2, grpnorm))))
      } else {
        norms <- apply(coefs, 2, grpnorm)
      }
      if(ncol(norms) > 1L){
        rownames(norms) <- object$term.labels
      } else {
        norms <- c(norms)
        names(norms) <- object$term.labels
      }
      return(norms)
    }
    
    
    
    ######***######   INITIAL CHECKS (CONTINUED)   ######***######
    
    ### check newx / newdata
    if(is.null(object$formula)){
      newx <- as.matrix(newx)
      ncoefs <- if(is.list(object$beta)) nrow(object$beta[[1]]) else nrow(object$beta)
      if(ncol(newx) != ncoefs) stop("Input 'newx' must satisfy:  ncol(newx) == nrow(object$beta)")
    } else {
      vnames <- rownames(attr(terms(object$formula), "factors"))[-1]
      check <- match(vnames, names(newdata))
      if(any(is.na(check))) stop("Input 'data' is missing variables included in model formula")
      if(!is.null(object$rk.args)){
        newx <- rk.model.matrix(object = object$formula, data = newdata,
                                knots = object$rk.args$knots,
                                Boundary.knots = lapply(object$rk.args$knots, function(x) x[c(1, length(x))]),
                                m = object$rk.args$m, 
                                periodic = object$rk.args$periodic,
                                xlev = object$rk.args$xlev)
      } else {
        newx <- model.matrix(object = object$formula, data = newdata)
      }
      if(colnames(newx)[1] == "(Intercept)") newx <- newx[,-1]
      ncoefs <- if(is.list(object$beta)) nrow(as.matrix(object$beta[[1]])) else nrow(as.matrix(object$beta))
      if(ncol(newx) != ncoefs) stop("Input 'newdata' produced a design matrix of the wrong dimension\n(likely due to a factor level mismatch between 'data' and 'newdata')")
      object$formula <- NULL
    }
    nobs <- nrow(newx)
    newxnames <- rownames(newx)
    if(is.null(newxnames)) newxnames <- 1:nobs
    
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
    
    
    
    ######***######   TERMS   ######***######
    
    if(type == "terms"){
      
      ### get coefs
      coefs <- coef(object, s = s)
      
      ### group names
      gnames <- object$term.labels
      ngroups <- object$ngroups
      group <- object$group
      
      ### remove intercept?
      int <- 0L
      if(gnames[1] == "(Intercept)") {
        int <- 1L
        gnames <- gnames[-1]
        ngroups <- ngroups - 1L
        group <- group[-1]
      }
      
      ### multinomial or other?
      if(family == "multinomial"){
        
        ## get number of response classes
        nresp <- length(object$beta)
        
        ## initialize terms list (of arrays)
        terms <- vector("list", nresp)
        names(terms) <- object$ylev
        for(j in 1:nresp){
          terms[[j]] <- array(dim = c(nobs, ngroups, ns))
          dimnames(terms[[j]]) <- list(newxnames, gnames, paste0("s", 1:ns))
        }
        
        ## loop thru terms
        for(k in 1:ngroups){
          id <- which(group == k)
          for(j in 1:nresp){
            terms[[j]][,k,] <- newx[,id,drop=FALSE] %*% coefs[[j]][id+int,,drop=FALSE]
          }
        }
        
        ## return results
        return(lapply(terms, drop))
        
      } else {
        
        ## initialize terms array
        terms <- array(dim = c(nobs, ngroups, ns))
        dimnames(terms) <- list(newxnames, gnames, paste0("s", 1:ns))
        
        ## loop thru terms
        for(k in 1:ngroups){
          id <- which(group == k)
          terms[,k,] <- newx[,id,drop=FALSE] %*% coefs[id+int,,drop=FALSE]
        } # end for(k in 1:object$ngroups)
        
        ## return results
        return(drop(terms))
        
      } # end if(family == "multinomial")
      
    } # end if(type == "terms")
    
    
    
    ######***######   IMPORTANCE   ######***######
    
    if(type == "importance"){
      
      ### get terms
      terms <- predict.grpnet(object, newx = newx, s = s, type = "terms")
      
      ### multinomial or other?
      if(family == "multinomial"){
        
        ## number of terms
        nterms <- dim(terms[[1]])[2]
        termnames <- dimnames(terms[[1]])[[2]]
        
        ## get number of response classes
        nresp <- length(object$beta)
        
        ## list (of matrices) for results
        imp <- vector("list", nresp)
        names(imp) <- object$ylev
        mat <- matrix(0, nrow = nterms, ncol = ns)
        rownames(mat) <- termnames
        colnames(mat) <- paste0("s", 1:ns)
        
        ## loop thru 1:nresp
        for(j in 1:nresp){
          
          imp[[j]] <- mat
          
          if(ns == 1L){
            etai <- scale(terms[[j]], scale = FALSE)
            eta0 <- rowSums(etai)
            imp[[j]] <- colSums(etai * eta0) / sum(eta0^2)
          } else {
            for(i in 1:ns){
              etai <- scale(terms[[j]][,,i], scale = FALSE)
              eta0 <- rowSums(etai)
              imp[[j]][,i] <- colSums(etai * eta0) / sum(eta0^2)
            } 
          } # end if(ns == 1L)
          
        } # end for(j in 1:nresp)
        
        if(ns == 1L) {
          imp <- sapply(imp, c)
          colnames(imp) <- object$ylev
        }
        return(imp)
        
      } else {
        
        ## number of terms
        nterms <- dim(terms)[2]
        termnames <- dimnames(terms)[[2]]
        
        ### matrix for results
        imp <- matrix(0, nrow = nterms, ncol = ns)
        rownames(imp) <- termnames
        colnames(imp) <- paste0("s", 1:ns)
        
        ## loop thru 1:ns
        if(ns == 1L){
          etai <- scale(terms, scale = FALSE)
          eta0 <- rowSums(etai)
          imp <- colSums(etai * eta0) / sum(eta0^2)
        } else {
          for(i in 1:ns){
            etai <- scale(terms[,,i], scale = FALSE)
            eta0 <- rowSums(etai)
            imp[,i] <- colSums(etai * eta0) / sum(eta0^2)
          } 
        } # end if(ns == 1L)
        
        return(drop(imp))
        
      } # end if(family == "multinomial")
      
    } # end if(type == "importance")
    
    
    
    ######***######   MULTINOMIAL   ######***######
    
    if(family == "multinomial"){ 
      
      # note: needs special code to return predictions for each response class
      
      ### get number of response classes
      nresp <- length(object$beta)
      
      ### initialize array for predictions
      fit <- array(data = NA, dim = c(nobs, nresp, ns))
      dimnames(fit) <- list(newxnames, object$ylev, paste0("s", 1:ns))
      
      ### predictions at object$lambda are easy...
      if(!newlambdas){
        for(k in 1:nresp) {
          if(nlam == 1L){
            fit[,k,] <- cbind(1, newx) %*% c(object$a0[k], object$beta[[k]])
          } else {
            fit[,k,] <- cbind(1, newx) %*% rbind(object$a0[k,], object$beta[[k]])
          }
        }
        if(type %in% c("response", "class")) {
          for(i in 1:ns) fit[,,i] <- object$family$linkinv(fit[,,i])
          if(type == "class") fit <- matrix(object$ylev[apply(fit, c(1,3), which.max)], nrow = nobs, ncol = ns)
        }
        return(drop(fit))
      }
      
      ### only 1 fit lambda?
      if(nlam == 1L){
        for(k in 1:nresp){
          fit[,k,] <- cbind(1, newx) %*% c(object$a0[k], object$beta[[k]])
        }
        if(type == "response") {
          fit <- object$family$linkinv(fit)
        } else if(type == "class"){
          # note: family == "multinomial" is implied
          fit <- object$ylev[apply(object$family$linkinv(fit), 1, which.max)]
        }
        return(fit)
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
      for(k in 1:nresp){
        coefs <- rbind(object$a0[k,], object$beta[[k]])
        coefs <- coefs[, left, drop=FALSE] %*% diag(sfrac, nrow = ns, ncol = ns) + coefs[, right, drop=FALSE] %*% diag(1 - sfrac, nrow = ns, ncol = ns)
        colnames(coefs) <- paste0("s", 1:ns)
        fit[,k,] <- cbind(1, newx) %*% coefs
      }
      
      ### return fitted values
      if(type %in% c("response", "class")) {
        for(i in 1:ns) fit[,,i] <- object$family$linkinv(fit[,,i])
        if(type == "class") fit <- matrix(object$ylev[apply(fit, c(1,3), which.max)], nrow = nobs, ncol = ns)
      }
      return(drop(fit))
      
    } # end if(object$family$family == "multinomial")
    
    
    
    ######***######   OTHER FAMILIES   ######***######
    
    ### predictions at object$lambda are easy...
    if(!newlambdas){
      if(nlam == 1L){
        fit <- cbind(1, newx) %*% c(object$a0, object$beta)
      } else {
        fit <- cbind(1, newx) %*% rbind(object$a0, object$beta)
      }
      if(type == "response") {
        fit <- object$family$linkinv(fit)
      } else if(type == "class"){ 
        # note: family == "binomial" is implied
        fit <- ifelse(object$family$linkinv(fit) <= 0.5, object$ylev[1], object$ylev[2])
      }
      return(drop(fit))
    }
    
    ### only 1 fit lambda?
    if(nlam == 1L){
      fit <- cbind(1, newx) %*% c(object$a0, object$beta)
      if(type == "response") {
        fit <- object$family$linkinv(fit)
      } else if(type == "class"){
        # note: family == "binomial" is implied
        fit <- ifelse(object$family$linkinv(fit) <= 0.5, object$ylev[1], object$ylev[2])
      }
      return(drop(matrix(fit, nrow = nrow(newx), ncol = ns)))
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
    colnames(coefs) <- paste0("s", 1:ns)
    
    ### return fitted values
    fit <- cbind(1, newx) %*% coefs
    if(type == "response") {
      fit <- object$family$linkinv(fit)
    } else if(type == "class"){
      # note: family == "binomial" is implied
      fit <- ifelse(object$family$linkinv(fit) <= 0.5, object$ylev[1], object$ylev[2])
    }
    return(drop(fit))
    
  } # end predict.grpnet