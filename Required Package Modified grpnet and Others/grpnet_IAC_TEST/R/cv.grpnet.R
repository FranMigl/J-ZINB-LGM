# Some code (for printing and plotting) is re-purposed from the R package
# "glmnet" (Hastie et al., 2010) https://cran.r-project.org/package=glmnet

cv.grpnet <-
  function(x, ...){
    UseMethod("cv.grpnet")
  } # end cv.grpnet 

print.cv.grpnet <- 
  function(x, digits = max(3, getOption("digits") - 3), ...){
    # convert to Title Case
    titlecase <- function(x){
      cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
      words <- strsplit(x, split = " ", fixed = TRUE)
      words <- sapply(words, cap1)
      paste(words, collapse = " ")
    }
    # enhanced measure label
    if(x$type.measure == "deviance"){
      if(x$grpnet.fit$family$family == "negative.binomial"){
        ylab <- paste("Negative Binomial", "Deviance")
      } else if(x$grpnet.fit$family$family == "inverse.gaussian"){
        ylab <- paste("Inverse Gaussian", "Deviance")
      } else {
        ylab <- paste(titlecase(x$grpnet.fit$family$family), "Deviance")
      }
    } else if(x$type.measure == "mse"){
      ylab <- "Mean Squared Error"
    } else if(x$type.measure == "mae"){
      ylab <- "Mean Absolute Error"
    } else if(x$type.measure == "class"){
      ylab <- "Misclassification Error"
    }
    # print call
    cat("\nCall:   ")
    print(x$call)
    # print measure
    cat("\nMeasure:   ", ylab, "\n\n")
    # print fit info
    df <- data.frame(Alpha = x$grpnet.fit$alpha,
                     Lambda = x$lambda[x$index], 
                     Index = x$index,
                     Measure = x$cvm[x$index],
                     SE = x$cvsd[x$index],
                     nzGroup = x$nzero[x$index],
                     nzCoef = x$grpnet.fit$nzcoef[x$index],
                     Df = x$grpnet.fit$df[x$index])
    rownames(df) <- c("min ", "1se ")
    print(df, digits = digits)
  } # end print.cv.grpnet

plot.cv.grpnet <-
  function(x, sign.lambda = 1, nzero = TRUE, ...){
    
    # check x
    if(!inherits(x, "cv.grpnet")) stop("Input 'x' must be an object of class 'cv.grpnet'")
    
    # check sign.lambda
    sign.lambda <- as.numeric(sign.lambda[1])
    if(!any(sign.lambda == c(-1, 1))) stop("Input 'sign.lambda' must be -1 or 1")
    xsign <- ifelse(sign.lambda < 0, "-", "")
    
    # define x coords and limits
    xcoords <- sign.lambda * log(x$lambda)
    xlimits <- range(xcoords)
    
    # enhanced y-axis label
    titlecase <- function(x){
      cap1 <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
      words <- strsplit(x, split = " ", fixed = TRUE)
      words <- sapply(words, cap1)
      paste(words, collapse = " ")
    }
    if(x$type.measure == "deviance"){
      if(x$grpnet.fit$family$family == "negative.binomial"){
        ylab <- "Negative Binomial Deviance"
      } else if(x$grpnet.fit$family$family == "inverse.gaussian"){
        ylab <- "Inverse Gaussian Deviance"
      } else {
        ylab <- paste(titlecase(x$grpnet.fit$family$family), "Deviance")
      }
    } else if(x$type.measure == "mse"){
      ylab <- "Mean Squared Error"
    } else if(x$type.measure == "mae"){
      ylab <- "Mean Absolute Error"
    } else if(x$type.measure == "class"){
      ylab <- "Misclassification Error"
    }
    
    # create args
    args <- list(...)
    args$x <- xcoords
    args$y <- x$cvm
    args$type = "n"
    if(is.null(args$ylim)) args$ylim <- c(min(x$cvlo), max(x$cvup))
    if(is.null(args$pch)) args$pch <- 20
    if(is.null(args$col)) args$col <- "red"
    if(is.null(args$xlab)) args$xlab <- if(sign.lambda < 0) expression("-" * log(lambda)) else expression(log(lambda))
    if(is.null(args$ylab)) args$ylab <- ylab
    
    # set up plot
    do.call(plot, args)
    
    # add nzeros
    if(nzero) axis(side = 3, at = xcoords, labels = x$nzero, tick = FALSE, line = -1)
    
    # add error bars
    width <- 0.01 * diff(xlimits)
    segments(xcoords, x$cvlo, xcoords, x$cvup, col = "darkgray")
    segments(xcoords - width, x$cvup, xcoords + width, x$cvup, col = "darkgray")
    segments(xcoords - width, x$cvlo, xcoords + width, x$cvlo, col = "darkgray")
    
    # add points
    points(x = args$x, y = args$y, pch = args$pch, col = args$col)
    
    # add lambda.min and lambda.1se
    abline(v = sign.lambda * log(x$lambda.min), lty = 3)
    abline(v = sign.lambda * log(x$lambda.1se), lty = 3)
    
  } # end plot.cv.grpnet