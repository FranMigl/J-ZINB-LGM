# Some code (for arguments and outputs) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

cv.grpnet.default <-
  function(x, 
           y, 
           group,
           weights = NULL,
           offset = NULL,
           alpha = c(0.01, 0.25, 0.5, 0.75, 1),
           type.measure = NULL,
           nfolds = 10, 
           foldid = NULL,
           same.lambda = FALSE,
           parallel = FALSE, 
           cluster = NULL, 
           verbose = interactive(), ...){
    # k-fold cross-validation for grpnet (default)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # Updated: 2023-09-04
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    ### get call
    cv.grpnet.call <- match.call()
    
    ### check x and y
    x <- as.matrix(x)
    nobs <- nrow(x)
    nvars <- ncol(x)
    ny <- if(is.matrix(y)) nrow(y) else length(y)
    if(ny != nobs) stop("Inputs 'x' and 'y' must satisfy:\nnrow(x) == length(y)  or  nrow(x) == nrow(y)")
    
    ### x names
    xnames <- colnames(x)
    if(is.null(xnames)) xnames <- paste0("x", 1:nvars)
    
    ### check group
    if(missing(group)) {
      group <- ingroup <- 1:nvars
      gsize <- rep(1L, nvars)
      names(gsize) <- 1:nvars
      ngrps <- nvars
    } else {
      if(length(group) != nvars) stop("Inputs 'x' and 'group' must satisfy:  ncol(x) == length(group)")
      ingroup <- group
      group <- as.factor(group)
      ngrps <- nlevels(group)
      if(ngrps > nvars) stop("Inputs 'x' and 'group' must satisfy:  ncol(x) >= nlevels(group)")
      gsize <- table(group)
      group <- as.integer(group)
    }
    
    ### check family
    args <- list(...)
    family <- args$family
    if(is.null(family)){
      family <- "gaussian"
    } else {
      family <- pmatch(as.character(family[1]), c("gaussian", "binomial", "multinomial", "poisson", "negative.binomial", "Gamma", "inverse.gaussian"))
      if(is.na(family)) stop("'family' not recognized")
      family <- c("gaussian", "binomial", "multinomial", "poisson", "negative.binomial", "Gamma", "inverse.gaussian")[family]
    }
    
    ### check weights
    if(is.null(weights)){
      noweights = TRUE
      weights <- rep(1.0, nobs)
    } else {
      noweights = FALSE
      weights <- as.numeric(weights)
      if(length(weights) != nobs) stop("Inputs 'y' and 'weights' must satisfy:  length(y) == length(weights)")
      if(any(weights < 0)) stop("Input 'weights' must be non-negative")
      weights <- weights / mean(weights)
    }
    
    ### check y
    ylev <- yfac <- NULL
    if(family == "gaussian"){
      y <- as.numeric(y)
    } else if(family == "binomial"){
      if(is.factor(y)){
        #if(nlevels(y) != 2L) stop("Input 'y' must be (coercible into) a factor with two levels when family = 'binomial'.")
        ylev <- levels(y)
        y <- ifelse(y == ylev[1], 0.0, 1.0)
      } else if(is.matrix(y)) {
        if(ncol(y) != 2L | any(y < 0)) stop("Input 'y' must be a matrix (# success, # failure) when family = 'binomial'")
        ytotal <- rowSums(y)
        weights <- weights * ytotal
        y <- y[,1] / ytotal
      } else {
        y <- as.numeric(y)
        if(any(y < 0) | any(y > 1)) stop("Input 'y' must contain values between 0 and 1 when family = 'binomial'" )
      }
      if(is.null(ylev)) ylev <- c(0, 1)
      yfac <- factor(ifelse(y <= 0.5, ylev[1], ylev[2]), levels = ylev)
    } else if(family == "multinomial"){
      if(is.factor(y)){
        yfac <- y
        yint <- as.integer(y)
        ylev <- levels(y)
        nlev <- nlevels(y)
        if(nlev < 3L) stop("Input 'y' must be (coercible into) a factor with 3 or more levels.")
        y <- matrix(0.0, nrow = nobs, ncol = nlev)
        colnames(y) <- ylev
        for(i in 1:nobs) y[i,yint[i]] <- 1.0
      } else if(is.matrix(y)){
        nlev <- ncol(y)
        ylev <- colnames(y)
        if(is.null(ylev)){
          ylev <- paste0("y", 1:nlev)
          colnames(y) <- ylev
        }
        y <- t(apply(y, 1, as.integer))
        if(min(y) < 0) stop("Input 'y' must be a matrix of counts (i.e, non-negative integers)")
        yrowsum <- rowSums(y)
        y <- y / yrowsum
        weights <- weights * yrowsum
        yfac <- factor(ylev[apply(y, 1, which.max)], levels = ylev)
      } else {
        stop("Invalid data input: 'y' must be a factor or matrix for when family = 'multinomial'.")
      }
    } else if (family == "poisson") {
      y <- as.integer(y)
      if(any(y < 0)) stop("Input 'y' must contain non-negative integers when family = 'poisson'.")
    } else if (family == "negative.binomial") {
      y <- as.integer(y)
      if(any(y < 0)) stop("Input 'y' must contain non-negative integers when family = 'negative.binomial'.")
    } else if(family == "Gamma") {
      y <- as.numeric(y)
      if(any(y <= 0)) stop("Input 'y' must contain non-negative numerics when family = 'Gamma'.")
    } else if(family == "inverse.gaussian") {
      y <- as.numeric(y)
      if(any(y <= 0)) stop("Input 'y' must contain non-negative numerics when family = 'inverse.gaussian'.")
    } # end if(family == "gaussian")
    
    ### check offset
    if(is.null(offset)){
      include.offset <- FALSE
      offset <- rep(0.0, nobs)
      if(family == "multinomial") offset <- matrix(offset, nrow = nobs, ncol = nlev)
    } else {
      include.offset <- TRUE
      if(family == "multinomial"){
        offset <- as.matrix(offset)
        if(nrow(offset) != nobs) stop("Inputs 'x' and 'offset' must satisfy:  nrow(x) == nrow(offset)")
        if(ncol(offset) == 1L){
          offset <- matrix(offset, nrow = nobs, ncol = nlev)
        } else {
          if(ncol(offset) != nlev) stop("Input 'offset' must be a vector of length nobs\n or a matrix of dimension nobs x nlev")
        }
      } else {
        offset <- as.numeric(offset)
        if(length(offset) != nobs) stop("Inputs 'y' and 'offset' must satisfy:  length(y) == length(offset)")
      }
    }
    
    ### check alpha
    alpha <- sort(unique(as.numeric(alpha)))
    if(any(alpha < 0 | alpha > 1)) stop("Input 'alpha' must satisfy:  0 <= alpha <= 1")
    nalpha <- length(alpha)
    
    ### check type.measure
    if(is.null(type.measure)){
      type.measure <- ifelse(family %in% c("binomial", "multinomial"), "class", "deviance")
    } else {
      type.measure <- pmatch(as.character(type.measure[1]), c("deviance", "mse", "mae", "class"))
      if(is.na(type.measure)) stop("Invalid 'type.measure' argument.")
      type.measure <- c("deviance", "mse", "mae", "class")[type.measure]
      if(type.measure == "class" && !(family %in% c("binomial", "multinomial")))
        stop("Input 'type.measure' can only be 'class' for binomial and multinomial family")
    }
    
    ### check nfolds and foldid
    if(is.null(foldid)){
      nfolds <- as.integer(nfolds[1])
      if(nfolds < 2L | nfolds > nobs) stop("Input 'nfolds' must satisfy:  2 <= nfolds <= nrow(x)")
      foldid <- sample(rep(1:nfolds, length.out = nobs))
    } else {
      foldid <- as.integer(foldid)
      if(length(foldid) != nobs) stop("Input 'foldid' must satisfy:  length(y) == length(foldid)")
      nfolds <- length(unique(foldid))
      if(nfolds < 2L | nfolds > nobs) stop("Input 'nfolds' must satisfy:  2 <= nfolds <= nrow(x)")
    } # end if(is.null(foldid))
    
    ### check parallel
    parallel <- as.logical(parallel[1])
    if(!any(parallel == c(TRUE, FALSE))) stop("Input 'parallel' must be TRUE or FALSE")
    
    ### check cluster
    if(parallel){
      if(nalpha == 1L) verbose <- FALSE
      if(is.null(cluster)){
        madeCluster <- TRUE
        #ncores <- parallel::detectCores()
        ncores <- 2L   # to ensure users don't accidentally spawn >2 processes
        cluster <- parallel::makeCluster(ncores)
      } else {
        madeCluster <- FALSE
        if(!any(class(cluster) == "cluster")) stop("Input 'cluster' must of an object of class 'cluster'")
      }
    } # end if(parallel)
    
    ### check verbose
    verbose <- as.logical(verbose[1])
    if(!any(verbose == c(TRUE, FALSE))) stop("Input 'verbose' must be TRUE or FALSE")
    
    
    ######***######   TUNE ALPHA   ######***######
    
    if(nalpha > 1L){
      alpha <- rev(alpha)
      alpha.tune <- data.frame(alpha = alpha, min.cvm = NA)
      # fit with alpha[1]
      if(verbose){
        cat("alpha = ", formatC(alpha[1], format = "f", digits = 2), ":", sep = "")
      }
      best <- cv.grpnet.default(x = x, 
                                y = y,
                                group = group,
                                weights = weights,
                                offset = offset,
                                alpha = alpha[1],
                                type.measure = type.measure,
                                nfolds = nfolds, 
                                foldid = foldid,
                                same.lambda = same.lambda,
                                parallel = parallel, 
                                cluster = cluster, 
                                verbose = FALSE,
                                ...)
      bestmin <- best$cvm[best$index[1]]
      alpha.tune[1,2] <- bestmin
      if(verbose){
        cat("   min(cvm) = ", bestmin, "\n", sep = "")
      }
      # loop thru alpha[k]
      for(k in 2:nalpha){
        if(verbose){
          cat("alpha = ", formatC(alpha[k], format = "f", digits = 2), ":", sep = "")
        }
        res <- cv.grpnet.default(x = x, 
                                 y = y,
                                 group = group,
                                 weights = weights,
                                 offset = offset,
                                 alpha = alpha[k],
                                 type.measure = type.measure,
                                 nfolds = nfolds, 
                                 foldid = foldid,
                                 same.lambda = same.lambda,
                                 parallel = parallel, 
                                 cluster = cluster, 
                                 verbose = FALSE,
                                 ...)
        resmin <- res$cvm[res$index[1]]
        alpha.tune[k,2] <- resmin
        if(verbose){
          cat("   min(cvm) = ", resmin, "\n", sep = "")
        }
        if(bestmin > resmin) {
          best <- res
          bestmin <- resmin
        }
      } # for(k in 2:nalpha)
      best$alpha <- alpha.tune
      best$call <- cv.grpnet.call
      best$grpnet.fit$ylev <- ylev
      return(best)
    } # end if(nalpha > 1L)
    
    
    ######***######   PRE-PROCESSING   ######***######
    
    ### collect fold ids
    fid <- vector("list", nfolds)
    for(k in 1:nfolds) fid[[k]] <- which(foldid == k)
    
    ### start progress bar
    if(verbose) pbar <- txtProgressBar(min = 0, max = nfolds + 1L, style = 3)
    
    ### generate lambda (if needed)
    grpnet.fit <- grpnet(x = x, 
                         y = y, 
                         group = group, 
                         weights = weights,
                         offset = offset, 
                         alpha = alpha,
                         ...)
    if(grpnet.fit$args$intercept){
      grpnet.fit$group <- c(0, ingroup)
    } else{
      grpnet.fit$group <- ingroup
    }
    grpnet.fit$ylev <- ylev
    lambda <- grpnet.fit$lambda
    nlambda <- length(lambda)
    if(verbose) setTxtProgressBar(pbar, 1)
    
    
    ######***######   K-FOLD CV   ######***######
    
    ### initialize matrix for results
    cvloss <- matrix(NA, nrow = nlambda, ncol = nfolds)
    
    ### separate work for multinomial family
    if(family == "multinomial"){
      
      if(parallel){
        
        # define parcvloss function
        parcvloss <-
          function(testid, xmat, ymat, group, family, weights, offset, alpha, nlambda,
                   lambda.min.ratio, lambda, penalty.factor, penalty, gamma, theta,
                   standardize, intercept, thresh, maxit, type.measure, same.lambda, yfac){
            temp <- grpnet(x = xmat[-testid,,drop=FALSE], 
                           y = ymat[-testid,,drop=FALSE], 
                           group = group, 
                           family = family,
                           weights = weights[-testid], 
                           offset = offset[-testid,],
                           alpha = alpha,
                           nlambda = nlambda,
                           lambda.min.ratio = lambda.min.ratio,
                           lambda = if(same.lambda) lambda else NULL,
                           penalty.factor = penalty.factor,
                           penalty = penalty,
                           gamma = gamma, 
                           theta = theta,
                           standardize = standardize,
                           intercept = intercept,
                           thresh = thresh,
                           maxit = maxit)
            temp$ylev <- levels(yfac)
            mu <- predict(temp, newx = xmat[testid,,drop=FALSE], 
                          s = if(same.lambda) NULL else lambda,
                          type = ifelse(type.measure == "class", "class", "response"))
            cvloss <- rep(NA, nlambda)
            if(type.measure == "deviance"){
              for(i in 1:nlambda) {
                cvloss[i] <- sum(temp$family$dev.resids(ymat[testid,,drop=FALSE], mu[,,i], weights[testid]))
              }
            } else if(type.measure == "mse") {
              for(i in 1:nlambda) {
                cvloss[i] <- mean((ymat[testid,,drop=FALSE] - mu[,,i])^2)
              }
            } else if(type.measure == "mae"){
              for(i in 1:nlambda) {
                cvloss[i] <- mean(abs(ymat[testid,,drop=FALSE] - mu[,,i]))
              }
            } else if(type.measure == "class"){
              cvloss <- 1 - colMeans(yfac[testid] == mu)
            }
            return(cvloss)
          } # end parcvloss
        
        # evaluate cvloss in parallel
        cvloss <- parallel::parSapply(cl = cluster, X = fid, FUN = parcvloss,
                                      xmat = x,
                                      ymat = y,
                                      group = group,
                                      family = family,
                                      weights = weights,
                                      offset = offset,
                                      alpha = grpnet.fit$alpha,
                                      nlambda = nlambda,
                                      lambda.min.ratio = lambda[nlambda] / lambda[1],
                                      lambda = lambda,
                                      penalty.factor = grpnet.fit$args$penalty.factor,
                                      penalty = grpnet.fit$args$penalty,
                                      gamma = grpnet.fit$args$gamma,
                                      theta = grpnet.fit$args$theta,
                                      standardize = grpnet.fit$args$standardize,
                                      intercept = grpnet.fit$args$intercept,
                                      thresh = grpnet.fit$args$thresh,
                                      maxit = grpnet.fit$args$maxit,
                                      type.measure = type.measure,
                                      same.lambda = same.lambda,
                                      yfac = yfac)
        
        # unvectorize
        cvloss <- matrix(cvloss, nrow = nlambda, ncol = nfolds)
        
      } else {
        
        for(k in 1:nfolds){
          if(same.lambda){
            temp <- grpnet(x = x[-fid[[k]],,drop=FALSE], 
                           y = y[-fid[[k]],,drop=FALSE], 
                           group = group, 
                           weights = weights[-fid[[k]]], 
                           offset = offset[-fid[[k]],,drop=FALSE],
                           alpha = alpha,
                           lambda = lambda, 
                           ...)
            temp$ylev <- ylev
            mu <- predict(temp, newx = x[fid[[k]],,drop=FALSE], 
                          type = ifelse(type.measure == "class", "class", "response"))
          } else {
            temp <- grpnet(x = x[-fid[[k]],,drop=FALSE], 
                           y = y[-fid[[k]],,drop=FALSE], 
                           group = group, 
                           weights = weights[-fid[[k]]], 
                           offset = offset[-fid[[k]],,drop=FALSE],
                           alpha = alpha,
                           ...)
            temp$ylev <- ylev
            mu <- predict(temp, newx = x[fid[[k]],,drop=FALSE], s = lambda, 
                          type = ifelse(type.measure == "class", "class", "response"))
          }
          if(type.measure == "deviance"){
            for(i in 1:nlambda) {
              cvloss[i,k] <- sum(grpnet.fit$family$dev.resids(y[fid[[k]],,drop=FALSE], mu[,,i], weights[fid[[k]]]))
            }
          } else if(type.measure == "mse") {
            for(i in 1:nlambda) {
              cvloss[i,k] <- mean((y[fid[[k]],,drop=FALSE] - mu[,,i])^2)
            }
          } else if(type.measure == "mae"){
            for(i in 1:nlambda) {
              cvloss[i,k] <- mean(abs(y[fid[[k]],,drop=FALSE] - mu[,,i]))
            }
          } else if(type.measure == "class"){
            cvloss[,k] <- 1 - colMeans(yfac[fid[[k]]] == mu)
          }
          if(verbose) setTxtProgressBar(pbar, k + 1)
        } # end for(k in 1:nfolds)
        if(verbose) close(pbar)
        
      } # end if(parallel)
      
    } else {
      
      if(parallel){
        
        # define parcvloss function
        parcvloss <- 
          function(testid, xmat, ymat, group, family, weights, offset, alpha, nlambda,
                   lambda.min.ratio, lambda, penalty.factor, penalty, gamma, theta,
                   standardize, intercept, thresh, maxit, type.measure, same.lambda, yfac){
            temp <- grpnet(x = xmat[-testid,,drop=FALSE], 
                           y = ymat[-testid], 
                           group = group, 
                           family = family,
                           weights = weights[-testid], 
                           offset = offset[-testid],
                           alpha = alpha,
                           nlambda = nlambda,
                           lambda.min.ratio = lambda.min.ratio,
                           lambda = if(same.lambda) lambda else NULL,
                           penalty.factor = penalty.factor,
                           penalty = penalty,
                           gamma = gamma,
                           theta = theta, 
                           standardize = standardize,
                           intercept = intercept,
                           thresh = thresh,
                           maxit = maxit)
            temp$ylev <- levels(yfac)
            mu <- predict(temp, newx = xmat[testid,,drop=FALSE], 
                          s = if(same.lambda) NULL else lambda,
                          type = ifelse(type.measure == "class", "class", "response"))
            if(type.measure == "deviance"){
              cvloss <- colSums(temp$family$dev.resids(ymat[testid], mu, weights[testid]))
            } else if(type.measure == "mse") {
              cvloss <- colMeans((ymat[testid] - mu)^2)
            } else if(type.measure == "mae"){
              cvloss <- colMeans(abs(ymat[testid] - mu))
            } else if(type.measure == "class"){
              cvloss <- 1 - colMeans(yfac[testid] == mu)
            }
            return(cvloss)
          } # end parcvloss
        
        # evaluate cvloss in parallel
        cvloss <- parallel::parSapply(cl = cluster, X = fid, FUN = parcvloss,
                                      xmat = x,
                                      ymat = y,
                                      group = group,
                                      family = family,
                                      weights = weights,
                                      offset = offset,
                                      alpha = grpnet.fit$alpha,
                                      nlambda = nlambda,
                                      lambda.min.ratio = lambda[nlambda] / lambda[1],
                                      lambda = lambda,
                                      penalty.factor = grpnet.fit$args$penalty.factor,
                                      penalty = grpnet.fit$args$penalty,
                                      gamma = grpnet.fit$args$gamma,
                                      theta = grpnet.fit$args$theta,
                                      standardize = grpnet.fit$args$standardize,
                                      intercept = grpnet.fit$args$intercept,
                                      thresh = grpnet.fit$args$thresh,
                                      maxit = grpnet.fit$args$maxit,
                                      type.measure = type.measure,
                                      same.lambda = same.lambda,
                                      yfac = yfac)
        
        # unvectorize
        cvloss <- matrix(cvloss, nrow = nlambda, ncol = nfolds)
        
      } else {
        for(k in 1:nfolds){
          if(same.lambda){
            temp <- grpnet(x = x[-fid[[k]],,drop=FALSE],
                           y = y[-fid[[k]]],
                           group = group,
                           weights = weights[-fid[[k]]],
                           offset = offset[-fid[[k]]],
                           alpha = alpha,
                           lambda = lambda, 
                           ...)
            temp$ylev <- ylev    # correct levels for classification loss
            mu <- predict(temp, newx = x[fid[[k]],,drop=FALSE],
                          type = ifelse(type.measure == "class", "class", "response"))
          } else {
            temp <- grpnet(x = x[-fid[[k]],,drop=FALSE],
                           y = y[-fid[[k]]],
                           group = group,
                           weights = weights[-fid[[k]]],
                           offset = offset[-fid[[k]]], 
                           alpha = alpha,
                           ...)
            temp$ylev <- ylev    # correct levels for classification loss
            mu <- predict(temp, newx = x[fid[[k]],,drop=FALSE], s = lambda,
                          type = ifelse(type.measure == "class", "class", "response"))
          } # end if(same.lambda)
          if(type.measure == "deviance"){
            cvloss[,k] <- colSums(grpnet.fit$family$dev.resids(y[fid[[k]]], mu, weights[fid[[k]]]))
          } else if(type.measure == "mse") {
            cvloss[,k] <- colMeans((y[fid[[k]]] - mu)^2)
          } else if(type.measure == "mae"){
            cvloss[,k] <- colMeans(abs(y[fid[[k]]] - mu))
          } else if(type.measure == "class"){
            cvloss[,k] <- 1 - colMeans(yfac[fid[[k]]] == mu)
          }
          if(verbose) setTxtProgressBar(pbar, k + 1)
        } # end for(k in 1:nfolds)
        if(verbose) close(pbar)
      } # end if(parallel)
      
    } # end if(family == "multinomial")
    
    
    ######***######   RESULTS   ######***######
    
    ### get mean and se
    cvm <- rowMeans(cvloss)
    cvsd <- apply(cvloss, 1, sd) / sqrt(nfolds)
    
    ### find lambda.min
    minid <- which.min(cvm)
    lambda.min <- lambda[minid]
    
    ### find lambda.1se
    min1se <- cvm[minid] + cvsd[minid]
    se1id <- which(cvm <= min1se)[1]
    lambda.1se <- lambda[se1id]
    
    ### close cluster
    if(parallel && madeCluster) parallel::stopCluster(cluster)
    
    ### correct term.labels
    if(grpnet.fit$args$intercept){
      grpnet.fit$term.labels <- c("(Intercept)", names(gsize))
    } else {
      grpnet.fit$term.labels <- names(gsize)
    }
    
    ### return results
    res <- list(lambda = lambda, 
                cvm = cvm, 
                cvsd = cvsd, 
                cvup = cvm + cvsd, 
                cvlo = cvm - cvsd,
                nzero = grpnet.fit$nzgrp, 
                grpnet.fit = grpnet.fit, 
                lambda.min = lambda.min,
                lambda.1se = lambda.1se,
                index = c(minid, se1id),
                type.measure = type.measure,
                call = cv.grpnet.call)
    class(res) <- "cv.grpnet"
    return(res)
    
  } # end cv.grpnet.default