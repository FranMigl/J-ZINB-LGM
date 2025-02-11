# Some code (for printing and plotting) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

grpnet <- 
  function(x, ...){
    UseMethod("grpnet")
  } # end grpnet

print.grpnet <-
  function(x, ...){
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    df <- data.frame(x$nzgrp,
                     x$df, 
                     round(100 * x$dev.ratio, 2),
                     x$lambda)
    colnames(df) <- c("nGrp", "Df", "%Dev", "Lambda")
    intercept <- ifelse(max(abs(x$a0)) > 0, TRUE, FALSE)
    print(df, digits = 6)
    cat("\n")
  } # end print.grpnet

plot.grpnet <-
  function(x, type = c("coef", "norm", "imp"), 
           newx, newdata, intercept = FALSE,
           color.by.group = TRUE, col = NULL, ...){
    int <- ifelse(x$args$intercept && !intercept, 1, 0)
    types <- c("coef", "norm", "imp")
    type <- pmatch(type[1], types)
    if(is.na(type)) stop("Invalid 'type' input")
    type <- types[type]
    if(color.by.group) {
      if(is.null(col)){
        col <- 1:x$ngroups
      } else {
        if(length(col) != x$ngroups) stop("Input 'col' must be of length x$ngroups")
      }
    } else {
      if(is.null(col)) col <- rep("black", x$ngroups)
    }
    if(type == "imp" && missing(newx) && missing(newdata)){
      stop("When type = 'imp', you need to provide either 'newx' or 'newdata'.")
    }
    res <- predict(x, newx = newx, newdata = newdata, type = type)
    if(type == "imp"){
      colors <- col
      if(x$family$family == "multinomial"){
        rnames <- rownames(res[[1]])
        cnames <- colnames(res[[1]])
        for(k in 1:length(res)){
          res[[k]] <- rbind(0, res[[k]])
          rownames(res[[k]]) <- c("(Intercept)", rnames)
          colnames(res[[k]]) <- cnames
        }
        index <- (1+int):nrow(res[[1]])
        for(j in 1:length(res)){
          plot(log(x$lambda), res[[j]][1,], ylim = extendrange(sapply(res, function(x) range(x[index,], na.rm = TRUE))),
               xlab = "Log Lambda", ylab = "Importance", t = "n", ...)
          legend("top", legend = x$ylev[j], bty = "n", cex = 0.8)
          for(k in index) {
            lines(log(x$lambda), res[[j]][k,], col = colors[k])
          }
        }
      } else {
        rnames <- rownames(res)
        cnames <- colnames(res)
        res <- rbind(0, res)
        rownames(res) <- c("(Intercept)", rnames)
        colnames(res) <- cnames
        index <- (1+int):nrow(res)
        plot(log(x$lambda), res[1,], ylim = extendrange(res[index,]),
             xlab = "Log Lambda", ylab = "Importance", t = "n", ...)
        for(k in index) {
          lines(log(x$lambda), res[k,], col = colors[k])
        }
      }
    } else if(type == "norm"){
      colors <- col
      index <- (1+int):nrow(res)
      plot(log(x$lambda), res[1,], ylim = extendrange(res[index,]),
           xlab = "Log Lambda", ylab = "L2 Norm", t = "n", ...)
      for(k in index) {
        lines(log(x$lambda), res[k,], col = colors[k])
      }
    } else {
      colors <- col[as.integer(as.factor(x$group))]
      if(x$family$family == "multinomial"){
        index <- (1+int):nrow(res[[1]])
        for(j in 1:length(res)){
          plot(log(x$lambda), res[[j]][1,], ylim = extendrange(sapply(res, function(x) range(x[index,], na.rm = TRUE))),
               xlab = "Log Lambda", ylab = "Coefficients", t = "n", ...)
          legend("top", legend = x$ylev[j], bty = "n", cex = 0.8)
          for(k in index) {
            lines(log(x$lambda), res[[j]][k,], col = colors[k])
          }
        }
      } else {
        index <- (1+int):nrow(res)
        plot(log(x$lambda), res[1,], ylim = extendrange(res[index,]),
             xlab = "Log Lambda", ylab = "Coefficients", t = "n", ...)
        for(k in index) {
          lines(log(x$lambda), res[k,], col = colors[k])
        }
      }
    }
    abline(h = 0)
  } # end plot.grpnet