# Here we have all the accessory functions to run the method

library(Matrix) # for sparse matrices

# The following function extract useful indexes to run the method
# In input requires only the list of k data matrices

indexes <- function(Xs, p){
  k = length(Xs)
  # Total number of rows and columns of the design matrix
  totrow <- sum(sapply(Xs, nrow)) 
  totcol <- sum(sapply(Xs, ncol))
  
  # Cumulative sums
  row_indices <- cumsum(sapply(Xs, nrow))
  col_indices <- cumsum(sapply(Xs, ncol))
  
  ##################################### Intercept matrix:
  
  pattern_matrix <- Matrix::Matrix(data = 0, nrow =  totrow, ncol = k, sparse = TRUE)
  for (j in 1:k) {
    pattern_matrix[row_indices[j] - nrow(Xs[[j]]) + 1:nrow(Xs[[j]]), j] <- 1
  }
  
  ##################################### Design Matrix without intercepts
  
  matricione <- Matrix(data = 0, nrow = totrow, ncol = totcol, sparse = TRUE)
  # nel "matricione" stiamo inserendo diagonalmente le k matrici di dati
  for (j in 1:k) {
    matricione[row_indices[j] - nrow(Xs[[j]]) + 1:nrow(Xs[[j]]),
               col_indices[j] - ncol(Xs[[j]]) + 1:ncol(Xs[[j]])] <- Xs[[j]]
  }
  
  
  # Final design matrix
  
  matricione <- cbind(pattern_matrix, matricione)
  
  
  nk <- sapply(Xs, nrow) # samples of each matrix
  ncsum = cumsum(nk) # cumulative number of samples
  inds <- Map(function(start, end) start:end, c(1, head(ncsum, -1) + 1), ncsum) 
  # k long list with to address row indexes per original matrix in the design matrix
  
  # grpnet requires group variables to be consecutives e.g. 111222333 etc.
    group = c(1:k, rep((k+1):(p+k), each = k)) 
  
  # Index that order the matrix according to the groups;
  
  ind <- as.numeric(c(1:k, sapply(1:(p), function(i){
    k + seq(from = i, to = k * p, by = p)
  })))
  
  # Same index withou intercepts
  ind_j <- as.numeric(sapply(1:(p), function(i){
    seq(from = i, to = k * p, by = p)
  }))
  
  # total number of column of the design matrix P
  
  P = ncol(matricione) # numero totale di colonne della design matrix
  
  output= list()
  output$matricione = matricione
  output$totrow = totrow
  output$totcol = totcol
  output$nk = nk
  output$ncsum = ncsum
  output$inds = inds
  output$group = group
  output$ind = ind
  output$ind_j = ind_j
  output$P = P
  return(output)
  # Useful indexes are:
  # totrow totcol: matrices indexes without intercepts
  # row_indices col_indices: dimensions increment for each matrix intercepts excluded
  # nk: non incremental matrices dimensions
  # ncsum: equal to row_indices (redundant)
  # inds: indexes to address each matrix in the design matrix
  # group: group indexes for grpnet
  # ind: to indicisize the matrix according to the group variable
  # ind_j: to fill coefficinets matrix excluding intercepts
  # P: final variables dimension of the design matrix
  # n, p = number of rows (genes), and variables (genes)
}

# ZILGM function to find lammax
find_lammax = function(X)
{
  tmp = t(X) %*% X
  lammax = 1/nrow(X) * max(abs(tmp[upper.tri(tmp)]))
  return(lammax)
}

# ZINB_GLGM function to find lammax
find_lammax_grp = function(X_list)
{
  # let's FAST compute T(X) %*% X
  Grammatrices <- lapply(X_list, crossprod)
  
  k = length(X_list)
  # Let's initialise the matrix containing the norms
  n = nrow(Grammatrices[[1]])
  
  norm_matrix <- matrix(0, nrow = n, ncol = n)
  
  # Let's compute the norms
  for (i in 1:n) {
    for (j in 1:n) {
      ijs <- sapply(Grammatrices, function(mat) mat[i, j])
      norm_matrix[i, j] <- sqrt(sum(ijs^2))
    }
  }
  # numsample <- norm(as.matrix(sapply(X_list, nrow)), type = "2")^2
  # numsample <- norm(as.matrix(sapply(X_list, nrow)), type = "2") # discuterne in riunione
  numsample <- sqrt(k)*sum(sapply(X_list, nrow)) # up to 23/09/24 the better
  # numsample <- mean(as.matrix(sapply(X_list, nrow)))
  rhomax = max(abs(norm_matrix[upper.tri(norm_matrix)]))
  lammax = 1/numsample * rhomax
  print(paste("Lammax: ", lammax, sep = ""))
  return(lammax)
}

# here we cast an alternative function related to finding the lammax
find_lammax_grp_alt = function(X_list)
{
  X <- do.call(rbind, args = X_list)
  tmp = t(X) %*% X
  rhomax = max(abs(tmp[upper.tri(tmp)]))
  print(paste("rhomax: ", rhomax, sep = ""))
  lammax = 1/nrow(X) * rhomax
  return(lammax)
}

# Here we have the function that transforms the grouped results coefficients in a joint network
per_lambda_normal <- function(output_matrix, indgrp, p, lambda) {
  # Initialize the result matrix for a single lambda
  norm_matrix_single <- matrix(0, nrow = p, ncol = p)
  
  # Loop over i and j to compute the norm for each element
  for (i in 1:p) {
    for (j in 1:p) {
      
      # Extract the corresponding coefficients from the output matrix for the current lambda
      coefficients <- output_matrix[indgrp[[i]],j,lambda]
      
      # Compute the norm of the extracted coefficients and store in the result matrix
      norm_matrix_single[i, j] <- sqrt(sum(coefficients^2))
    }
  }
  
  return(norm_matrix_single)
}


# Negative Binomial density function
dNBI = function(y, mu, theta, log = FALSE)
{
  density = lgamma(y + theta + 1e-10) - lgamma(y + 1) - lgamma(theta + 1e-10) + theta * (log(theta + 1e-10) - log(theta + mu + 1e-10)) + y * (log(mu + 1e-10) - log(theta + mu + 1e-10))
  if (log == FALSE) {density = exp(density)}
  return(density)
}

# NB second parametrisation; that's the first parametrisation where sigma = alpha = mu / theta.
#  Infatti theta = mu / sigma in dNBI comporta  mu * Theta / mu = Theta
dNBII = function(y, mu, sigma, log = FALSE) 
{
  density = dNBI(y, mu = mu, theta = mu / sigma, log = log)
  return(density)
}

# PSEUDOLIKELIHOOD for the second parametrisation
grp_nb2_objective = function(y, weights, prob, bvec, mu, sigma = NULL, lambda, penalty.factor, posz, k, p)
{
  # Firstly we divide the beta results in groups, intercepts excluded
  # These will be the values to be l2 norm computed for the penalty
  bgrp <- lapply(seq(k + 1, dim(bvec)[1], by = k), function(start) {
    matrix(bvec[start:(start + k - 1), ])
  })
  
  # Iteratively computing penalty
  penalty = 0
  for (i in 1:p){
    penalty = penalty + penalty.factor[i + k]*norm(bgrp[[i]])
  }
  
  # actual likelihood computation
  pnl = 0
  for (i in 1:k) {
    pnl = pnl - sum(weights[inds[[i]]] * log(prob[i] * posz[inds[[i]]] + (1 - prob[i]) * dNBII(y = y[inds[[i]]], sigma = sigma, mu = mu[inds[[i]]],log = FALSE) + 1e-10))
  }
  return(pnl/k + lambda*penalty)
}

# Theta estimate
theta_ml = function(y, mu, w = NULL) {
  n = length(y)
  if (is.null(w)) {w = rep(1, n)}
  nb_theta = function(theta, mu, y, w) {
    return(sum(w * dNBI(y = y, theta = theta, mu = mu, log = TRUE)))
  }
  fit = optimize(nb_theta, y = y, mu = mu, w = w, interval = c(1e-4, 5e+3), maximum = TRUE)
  theta = ifelse(fit$maximum > 1e+3, 1e+8, fit$maximum)
  return(theta)
}

# Sigma estimate, i.e. theta for NB2 parametrisation
sigma_ml = function(y, mu, weights = NULL)
{
  NB2_theta = function(sigma, mu, y, weights) {
    return(sum(weights * dNBII(y = y, sigma = sigma, mu = mu, log = TRUE)))
  }
  # start = c(0.01)
  sigma = tryCatch(
    {result = optimize(NB2_theta, y = y, mu = mu, weights = weights, interval = c(1e-6, 1000), maximum = TRUE)
    sigma = result$maximum
    return(sigma)
    },
    warning = function(w) {
      cat("Warning message:", conditionMessage(w), "\n")
      sigma = (mean(y)^2)/(var(y)-mean(y))
      print(paste("Sigma after warning :", sigma, sep = " "))
      return(sigma)
    }, error = function(e) {
      cat("Error message:", conditionMessage(e), "\n")
      sigma = (mean(y)^2)/(var(y)-mean(y))
      print(paste("Sigma after Error :", sigma, sep = " "))
      return(sigma)
    })
  
  # If an error or warning occurred, result will contain the value returned by the error or warning handler.
  # If no error or warning occurred, result will contain the result of the code block.
  
  # sigma = ifelse(sigma <= 5e-5, 0, sigma)
  return(sigma)
}

# Network computing function
hat_net = function(coef_mat, thresh = 1e-6, type = c("AND", "OR"))
{
  type = match.arg(type)
  
  tmp_mat = abs(coef_mat) > thresh
  
  if (type == "AND") {
    res_mat = tmp_mat * t(tmp_mat)
  }
  
  if (type == "OR") {
    res_mat = (tmp_mat + t(tmp_mat) > 0) * 1
  }
  return(res_mat)
}

# The following is a function to compute performances, once the model is run

metrics = function (net, hat) # Works on adjacency matrices hat, given a network net
{
  net <- as_adjacency_matrix(net, type = "both", sparse = FALSE) # diagonal values zero
  hat = abs(hat != 0) # binarize the zero diagonal matrix
  p = ncol(net)
  pfms = matrix(0, p, 6)
  colnames(pfms) = c("ZE","TN", "FP", "NZ", "FN", "TP")
  rownames(pfms) = 1:p
  
  for (j in 1:p) {
    # negative cases
    flag = net[, j] == 0 # logical vector of zeroes per column
    flag[j] = FALSE # Excluding the zero diagonal elements from the logical vector
    pfms[j, 1] = sum(flag)/2 # ZE zeroes elements in net 
    pfms[j, 2] = sum(net[flag, j] == hat[flag, j])/2 # TN
    pfms[j, 3] = sum(net[flag, j] != hat[flag, j])/2 # FP
    
    # positive cases
    flag = net[, j] != 0 # logical vector of non zero values
    pfms[j, 4] = sum(flag)/2 # NZ Non zeroes elements in net
    pfms[j, 5] = sum(net[flag, j] != hat[flag, j])/2 # FN
    pfms[j, 6] = sum(net[flag, j] == hat[flag, j])/2 # TP
    
  }
  
  x <- colSums(pfms)
  x[7] = x[6] / x[4]
  x[8] = x[5] / x[4]
  x[9] = x[2] / x[1]
  x[10] = x[3] / x[1]
  x[11] = x[6] /(x[6] + x[3]) 
  x[12] = (2*x[6]*x[2] + x[6]*x[3] + x[5]*x[2])/(2*(x[6]+x[5])*(x[2]+x[3]))
  x[13] = 2*(x[11] * x[7])/(x[11] + x[7])
  x[14] = round((x[3] + x[6])) 
  names(x)[7:14] <- c("TPR", "FNR", "TNR", "FPR", "Precision", "Balanced_Accuracy", "F1", "N.edges")
  x <- x[c(1:6,14,7:13)]
  return(x)
}