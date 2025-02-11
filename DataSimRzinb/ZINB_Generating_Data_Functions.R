# Function generating one matrix according to a given graph

Generetwork_NB <- function(net, # ground truth
                           seed, # to have the same matrix if desired
                           n, # how many samples
                           p, # how many variables (graph size)
                           theta, # overdispersion of sources
                           means, # NB sources mean
                           zerop # zero percentage
                           )
  {
  # We'll have the graph realisation X as Y%*%B, so let's start from the dimension of matrix Y; an n row times p + p(p-1)/2 columns sources matrix
  bino <- choose(p, 2) # p(p-1)/2
  m <- p + bino # p + p(p-1)/2  
  
  # We set the seed
  set.seed(seed)
  matseed <- .Random.seed
  
  # We can now generate Y matrix values:
  Y <- matrix(emdbook::rzinbinom(n = n * m,
                            mu = means, # means from park chosen 1.5
                            size = theta,
                            zprob = zerop),
              nrow = n,
              ncol = m) 
  
  E <- matrix(emdbook::rzinbinom(n = n*p, mu = 0, size = 1e8, zprob = zerop), nrow = n, ncol = p)
  
  # We now got to build matrix B, this will be made by an identity matrix p times p, a permutation matrix and the vectorization of the upper triangular matrix of the graph adjacency matrix
  
  Ip <- diag(p)
  
  # To generate the permutation matrix let's start from all the possible two element combination
  combinations <- combn(1:p, 2)
  
  
  # And then we can get the permutation matrix
  permutation_matrix <- matrix(0, nrow = p, ncol = choose(p, 2))
  for (i in 1:ncol(permutation_matrix)) {
    permutation_matrix[combinations[1,i], i] <- 1
    permutation_matrix[combinations[2,i], i] <- 1
  }
  
  # Finally we retrieve the upper triangular vectorization of the adjacency matrix by row order
  
  # We exploit the fact that under diagonal elements of a symmetric matrix taken by column correspond to the upper one take by row
  adj_matrix <- igraph::as_adjacency_matrix(net, type = "lower", sparse = FALSE)
  
  # With these elements we can now build the fundamental vector giving our network structure to the data
  
  upper_tri_vector <- adj_matrix[lower.tri(adj_matrix)]
  upper_tri_matrix <- t(replicate(p, upper_tri_vector))
  
  # We cast the graph structure into the data by taking the product between the permutation matrix and the upper triangular matrix
  
  K <-permutation_matrix * upper_tri_matrix
  
  # Finally, to get matrix b, we vertically stuck K under the identity matrix [Ip]=pxp
  
  B <- t(cbind(Ip, K))
  
  # Our X matrix, simulating a realisation of the graph, possibly representing a scRNA-seq cluster is:
  
  X <- Y%*%B + E
  
  return(X = X)
}

# Function to generate multiple matrices from the same network 

k_net_datagen <- function(net,
                          k,
                          seed,
                          n,
                          p,
                          theta,
                          means,
                          zerop)
  {
  grpnetworks <- list()
  for (i in 1:k){
    grpnetworks[[i]] <- Generetwork_NB(net = net, seed = seed + i,
                                       n, p, theta, means = means[i], zerop = zerop)
    # net indica quale rete usare per la struttura dei dati
    # tramite seed, fissata una rete, la matrice risultante cambia
    # n è la sample size della matrice generata, ogni vertice produce n osservazioni
    # p è il numero di colonne ed equivale al numero di vertici del grafo
    # Theta controlla la dispersione mediante la formula var = mean + theta*mean^2
    # means controlla le medie delle sorgenti, vedi View(Generetwork_NB_grpnet)
  }
  return(grpnetworks)
}

# Function

Dropout_Injection <- function(X,
                              dropout_rate,
                              source,
                              seed = 1926) 
  {
  X_dropout <- X
  set.seed(seed)
  X_dropout[X_dropout < 5*source] <- ifelse(runif(length(X_dropout[X_dropout < 5*source])) <= dropout_rate,
                                            0, X_dropout[X_dropout < 5*source])
  return(X_dropout)
}
Dropout_Injection_meandeg <- function(X,
                              dropout_rate,
                              source,
                              seed = 1926,
                              meandeg) 
{
  X_dropout <- X
  set.seed(seed)
  X_dropout[X_dropout < 5*source*meandeg] <- ifelse(runif(length(X_dropout[X_dropout < 5*source*meandeg])) <= dropout_rate,
                                            0, X_dropout[X_dropout < 5*source*meandeg])
  return(X_dropout)
}

# Dropout calculation
calculate_dropout <- function(X) {
  zeroes <- sum(X == 0) / length(X)
  return(zeroes)
}