# The first function just run the method over temporary lambdas, if it is not provided.
# These provisory lambdas then allow to obtain a grid over which the method  is run
# As output one get the networks and other stuff

zinb_LGM_grp <- function(
    Xlist,
    lambda = NULL,
    nlambda = 50,
    conperc = 50,
    sym = c("AND", "OR"),
    offset = NULL,
    theta = NULL,
    thresh = 1e-6,
    weights_mat = NULL,
    penalty_mat = NULL,
    nCores = 1,
    verbose = 0
) {
  # Match the symmetry argument
  sym <- match.arg(sym)
  fun_call <- match.call()
  
  # Initialize offset if null
  if (is.null(offset)) {
    offset <- unlist(lapply(Xlist, function(i) log(rowSums(i))))
  }
  
  nvar <- dim(Xlist[[1]])[2]
  if (nvar < 2) stop("X must be a matrix with 2 or more columns")
  
  # Create indexes for matrices
  mat_indexes <- indexes(Xlist, nvar)
  ratio <- ifelse(any(sapply(Xlist, function(i) dim(i)[1] < nvar)), 0.05, 0.0001)
  penalty <- "LASSO"
  
  # The function that computes the number of edge for each selected lambda in the bisection process
  adjust_lambda <- function(current_lambda) {
    tmp_net <- zinb_LGM_net_grp(
      X = mat_indexes$matricione,
      lambda = current_lambda,
      sym = sym,
      theta = theta,
      thresh = thresh,
      weights_mat = weights_mat,
      penalty_mat = penalty_mat,
      nCores = nCores,
      p = nvar,
      totrow = mat_indexes$totrow,
      totcol = mat_indexes$totcol,
      P = mat_indexes$P,
      k = length(Xlist),
      nk = mat_indexes$nk,
      inds = mat_indexes$inds,
      group = mat_indexes$group,
      ind = mat_indexes$ind,
      ind_j = mat_indexes$ind_j,
      offset = offset,
      verbose = 0
    )
    
    edges <- sum(unlist(lapply(tmp_net$hat_net, function(x) {
      igraph::gsize(igraph::graph_from_adjacency_matrix(as(x, "CsparseMatrix"), mode = "undirected"))
    })))
    list(edges = edges, tmp_net = tmp_net)
  }
  
  # Lambda computation
  if (is.null(lambda)) {
    if (verbose > 0) cat("\nSearching lambda...\n")
    maxcon <- nvar * (nvar - 1) / 2
    lambda_max <- find_lammax_grp(Xlist)
    
    # Adjust lambda_max
    repeat {
      res <- adjust_lambda(lambda_max)
      if (res$edges >= (maxcon / 100)) break
      lambda_max <- lambda_max / 2
    }
    lambda_min <- lambda_max * ratio
    # Adjust lambda_min
    repeat {
      res <- adjust_lambda(lambda_min)
      if (res$edges <= (maxcon / 100) * conperc) break #res$edges >= (maxcon / 100) * conperc
      lambda_min <- lambda_min * 2 # lambda_min <- lambda_min / 2
    }
    
    # Generate sequence of lambda
    lambda <- exp(seq(log(lambda_min), log(lambda_max), length.out = nlambda))
    if (verbose > 0) cat("Lambda search complete\n")
  } else {
    nlambda <- length(lambda)
  }
  
  # Main computation
  if (verbose > 0) {
    cat(
      "Learning for NBII graphical model\n",
      "nlambda: ", nlambda, "\n",
      "penalty function: ", penalty, "\n",
      sep = ""
    )
  }
  
  net <- zinb_LGM_net_grp(
    X = mat_indexes$matricione,
    lambda = lambda,
    sym = sym,
    theta = theta,
    thresh = thresh,
    weights_mat = weights_mat,
    penalty_mat = penalty_mat,
    nCores = nCores,
    p = nvar,
    totrow = mat_indexes$totrow,
    totcol = mat_indexes$totcol,
    P = mat_indexes$P,
    k = length(Xlist),
    nk = mat_indexes$nk,
    inds = mat_indexes$inds,
    group = mat_indexes$group,
    ind = mat_indexes$ind,
    ind_j = mat_indexes$ind_j,
    offset = offset,
    verbose = verbose
  )
  
  if (verbose > 0) cat("Full process complete\n\n")
  
  # Output results
  list(
    network = net$hat_net,
    coef_network = net$coef_net,
    lambda = lambda,
    call = fun_call,
    itermat = net$itermat,
    warnings = net$warnings
  )
}

# The second function run the local grouped graphical model edges estimation parallely over each
# variable "j"

zinb_LGM_net_grp = function(X,
                            lambda = NULL,
                            sym = c("AND", "OR"),
                            theta = NULL,
                            thresh = 1e-6,
                            weights_mat = NULL,
                            penalty_mat = NULL,
                            nCores = 1,
                            p,
                            totrow,
                            totcol,
                            P,
                            k,
                            nk,
                            inds,
                            group,
                            ind,
                            ind_j,
                            offset = NULL,
                            verbose = 0)
  {
  
  nlambda = length(lambda) # salviamo la dimensione della griglia
  
  # la seguente matrice coef_mat deve accomodare i coefficienti risultanti escluse le intercette quindi saranno da ciascun fit k*p, per ogni gene quindi per p, per il numero di lambda nlambda 
  coef_mat = array(dim = c(k*p, p, nlambda))
  
  # creiamo anche una lista ove possiamo salvare gli eventuali warning che si possono aver avuto dai fit con il corrispettivo valore di j
  
  warnings <- list()
  
  # Inizializziamo una matrice di p colonne che accomodi i risultati del numero di iterazioni ottenute per ogni lambda, questa ci serve come controllo per capire quante iterazioni fa il nostro metodo
  
  itermat <- list() # matrix(nrow = nlambda, ncol = p)
  
  # weights_mat, è la possibile matrice dei pesi sulle osservazioni. Volendo si può specificare, altirmenti ogni osservazione avrà il medesimo peso. 
  # NB non si devono moltiplicare anche le intercette quindi ne sono p*k
  if (is.null(weights_mat)) {
    weights_mat = matrix(1, totrow, p*k)
  }
  # se i pesi sono forniti esternamente facciamo dei check
  if (any(weights_mat < 0)) {"The elements in weights_mat must have non-negative values"}
  if ((NROW(weights_mat) != totrow) | (NCOL(weights_mat) != totcol)) {"The number of elements in weights_mat not equal to the number of rows and columns on X, intercepts aside"}
  
  # Ora lanciamo la parallelizzazione
  
  # inizializziamo i cores
  cl <- parallel::makeCluster(nCores)
  
  
  
  # esportiamo variabili necessarie al ciclo
  parallel::clusterExport(cl, c("k", "p", "ind",
                                "ind_j", "group", "nlambda",
                                "verbose", "thresh", "X",
                                "lambda", "theta", "nk",
                                "inds", "weights_mat", "penalty_mat",
                                "P", "zinb_LGM_wrapper_grp",
                                "dNBI", "dNBII",
                                "grp_nb2_objective",
                                "sigma_ml", "theta_ml",
                                "zilgm_negbin2_grp", "hat_net", "per_lambda_normal"), environment())
  
  # carichiamo librerie necessarie
  parallel::clusterEvalQ(cl, {
    library(Matrix)
    library(grpnet)
    library(igraph)
    library(MASS)
    library(parallel)
  })
  
  # passiamo esplicitamente gli argomenti a cluster apply
  coef_tmp = parallel::clusterApply(cl = cl, 1:p, function(j,
                                                           group,
                                                           X,
                                                           lambda,
                                                           offset,
                                                           theta,
                                                           nk,
                                                           inds,
                                                           weights_mat,
                                                           penalty_mat,
                                                           k,
                                                           P,
                                                           ind,
                                                           ind_j,
                                                           p,
                                                           nlambda,
                                                           thresh,
                                                           verbose) {
    # Prendiamo gli indici di Y
    jth = j + seq(k, k*p, by = p) # indice j di ogni matrice nel matricione
    ind_1 = ind[-which(ind == jth)] # indice di tutti tranne i jth su cui fare regressione
    group_j = group[-which(group == j+k)] # stesso discorso togliamo i jth dai gruppi
    # lanciamo la seconda shell del nostro metodo: zilgm_wrapper_grp che a j fissato cicla su lambda. Questa funzione contiene anche l'implementazione di active set
    
    # NB penalty.factor è NULL, lo useremo per introdurre il termine della penalty di gruppo sqrt(k) più avanti ma va considerato variabile interna
    zinb_LGM_wrapper_grp(j = j,
                         X = X,
                         jth = jth,
                         lambda = lambda,
                         theta = theta,
                         nk = nk,
                         inds = inds,
                         group = group_j,
                         weights = weights_mat[, j],
                         penalty.factor = penalty_mat[, j],
                         k = k,
                         P = P,
                         ind = ind_1,
                         ind_j = ind_j,
                         offset = offset,
                         p = p,
                         nlambda = nlambda,
                         thresh = thresh,
                         verbose = verbose)
  },          group = group,
  X = X,
  lambda = lambda,
  theta = theta,
  nk = nk,
  inds = inds,
  weights_mat = weights_mat,
  penalty_mat = penalty_mat,
  k = k,
  P = P,
  ind = ind,
  ind_j = ind_j,
  p = p,
  nlambda = nlambda,
  offset = offset,
  thresh = thresh,
  verbose = verbose)
  
  
  parallel::stopCluster(cl)
  
  # Qui andiamo a popolare la matrice dei coefficienti e quelle dei warnings e iterazioni
  for (j in 1:p) {
    itermat[[j]] <- coef_tmp[[j]]$niterations
    coef_mat[, j, ] = as.matrix(coef_tmp[[j]]$Bmat)
    warnings[[j]] <- j
    if(!is.null(coef_tmp[[j]]$warning_lambda)){
      # Sostituisci la tua funzione di generazione di vettore casuale qui
      warnings[[j]] <- c(warnings[[j]],coef_tmp[[j]]$warning_lambda)
    }
  }
  
  # NB coef_mat non è come in Park, direttamente la matrice (p)x(p), è (k*p)x(p) e va trasformato
  # usiamo questo indice di gruppo
  indgrp <- split(ind_j, group[-(1:k)])
  
  # Conserveremo i risultati nella seguente matrice:
  coef_fin <- array(dim = c(p, p, nlambda))
  
  # riempiamo la matrice coef_fin
  for (l in 1:nlambda) {
    # Compute and store the matrix for the current lambda
    coef_fin[,,l] <- per_lambda_normal(coef_mat, indgrp, p, l)
  }
  
  
  # Futuro codice per costruire le reti:
  ghat = lapply(1:nlambda, FUN = function(l) hat_net(coef_fin[, , l], thresh = thresh, type = sym))
  gs = lapply(1:nlambda, FUN = function(l) Matrix(ghat[[l]]))
  
  # componiamo la lista risultato
  result_list <- list(coef_net = coef_mat, hat_net = gs, warnings = warnings, itermat = itermat)
  #return(list(hat_net = gs, coef_net = coef_mat)) # per quando avremo direttamente la rete in uscita, per ora i coefficienti
  return(result_list)
}

# The third function retrieve the local GM edges of a variable j for each lambda in the grid

zinb_LGM_wrapper_grp = function(j, X,
                                jth,
                                lambda, theta,
                                nk, inds, group,
                                weights, penalty.factor,
                                k, P, ind, ind_j, #aggiunto
                                offset,
                                n, p, nlambda, thresh,
                                verbose)
  {
  # j è FISSATO, CICLIAMO SU LAMBDA
  
  # la seguente matrice Bmat deve accogliere i coefficienti del fit ad ogni iterazione j e per ogni lambda.
  
  # la matrice Bmat quindi serve per accogliere i coefficienti ottenuti dal fit, ma anche quelli nulli diagonali che non si ottengono dal fit, i quali devono rimanere infatti nulli
  # b0 invece accoglierà le intercette
  
  # e.g. per p = 500, dal fit, quindi escludendo le intercette, ne otteniamo 1497 a gruppi i.e. 123 = Interc, 456 = 1° gruppo, 789 2° gruppo etc...
  
  # Bmat ne deve sempre accomodare 1500
  
  # Bmat = Matrix(0, p, nlambda, sparse = TRUE)
  Bmat = Matrix::Matrix(0, P - k, nlambda, sparse = TRUE) # Matrix matrice sparsa
  
  # b0 deve sempre accomodare i primi k (le intercette)
  b0 = Matrix::Matrix(0, k, nlambda, sparse = TRUE)
  
  # usiamo il seguente vettore per salvare il beta stimato dall'EM a ciascun lambda e passarlo come beta iniziale nel ciclo ai lambda successivi
  
  betaprev = NULL
  
  # Inizializziamo il vettore atto a raccogliere i warning
  warning_lambda <- c()
  
  # Inizializziamo anche il vettore preposto a raccogliere il numero di iterazioni per un dato lambda, fissato j
  niterations <- c()
  
  if (length(ind) == 0) { # Se non ci sono colonne selezionate restituisci risultato vuoto
    Bmat = Bmat
    b0 = b0
  } else {
    if(k == 1){
      y = X[,jth]
    } else {
      y = rowSums(X[,jth])
    }
    for (iter in 1:nlambda) {
      # Primo lambda: calcoliamo tutto agnosticamente
      if (iter == 1){
        # printiamo sempre if (verbose == 1) {
        cat("lambda = ", lambda[iter], ", ", j, "/", p, "th node learning \n", sep = "")
        # }
        
        # tryCatch per gestire i lambda problematici, ossia dove eventualmente il fit da errori
        coef_res <- tryCatch({
          fitem = zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE],
                                    k, 
                                    nk, inds, group,
                                    lambda = lambda[iter], theta = theta,
                                    bvec0 = NULL,
                                    w = weights,
                                    penalty.factor = penalty.factor,
                                    offset = offset,
                                    thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
        }, warning = function(w) {
          # salviamo quale lambda sta eventualmente dando problemi per questo j
          cat("Warning NAN nel fit su lambda:", lambda[iter], "\n")
          
          # Restituiamo comunque il risultato del fit
          fitem <- zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE],
                                     k, 
                                     nk, inds, group,
                                     lambda = lambda[iter], theta = theta,
                                     bvec0 = NULL,
                                     w = weights,
                                     penalty.factor = penalty.factor,
                                     offset = offset,
                                     thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
          fitem$warning_lambda = lambda[iter]
          return(fitem)
        })
        if(!is.null(coef_res$warning_lambda)){
          warning_lambda = c(warning_lambda, coef_res$warning_lambda)
        }
        
        # riempiamo solo con i coefficienti calcolati, lasciamo nulli gli altri
        Bmat[ind_j[ind_j != (jth-k)],iter] = coef_res$bvec[-(1:k)]
        
        # Aggiungiamo una condizione di stop
        if (sum(Bmat[,iter] == 0) == length(Bmat[,iter])){
          break}
        
        # Gli active sono proprio quelli diversi da zero in Bmat,dobbiamo però "tradurli" nell'ordine del fit di y VS X
        # Non usiamo Active, o almeno manteniamo la taglia del sistema
        # Active <- as.numeric(c(1:k, sapply(which(Bmat[1:p, iter] != 0), function(i){
          #k + seq(from = i, to = k * p, by = p)
        #})))
        
        # Non aggiorniamo dunque i gruppi del prossimo fit: in quanto vogliamo
        #  conservare la taglia del problema
        # group = c(1:k, rep((k+1):(k+length(which(Bmat[1:p, iter] != 0))), each = k))
        
        # Salviamo il beta nuovo da dare in input al prossimo step:
        # NO! betaprev = c(coef_res$bvec[1:k], coef_res$bvec[-(1:k)][coef_res$bvec[-(1:k)]!=0])
        # In questo caso vogliamo conservare lo stesso identico beta ottenuto dal fit 
        
        # Ulteriore modifica : non salviamo alcun betaprev, vogliamo che sia sempre vuoto
        # la dimensionalità deve restare invariata
        # betaprev = coef_res$bvec 
        
        # Riempiamo anche il vettore delle intercette 
        b0[,iter] = coef_res$bvec[1:k]
        
        # E salviamo quante iterazioni sono state necessarie
        niterations[iter] <- coef_res$iterations
        
      } else { # le altre iterazioni partono con un beta noto
        # printiamo sempre if (verbose == 1) {
        cat("lambda = ", lambda[iter], ", ", j, "/", p, "th node learning \n", sep = "")
        # }
        
        # tryCatch per gestire i lambda problematici
        coef_res <- tryCatch({
          fitem = zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE], # ind è cambiato con gli active
                                    k, 
                                    nk, inds, group,
                                    lambda = lambda[iter], theta = theta,
                                    bvec0 = NULL,  
                                    w = weights,
                                    penalty.factor = penalty.factor,
                                    offset = offset,
                                    thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
        }, warning = function(w) {
          # salviamo quale lambda sta dando problemi per questo j
          cat("Warning NAN nel fit su lambda:", lambda[iter], "\n")
          
          # Restituiamo comunque fitem
          fitem <- zilgm_negbin2_grp(y = y, x = X[, ind, drop = FALSE], # ind è cambiato con gli active
                                     k, 
                                     nk, inds, group,
                                     lambda = lambda[iter], theta = theta,
                                     bvec0 = NULL,
                                     w = weights,
                                     penalty.factor = penalty.factor,
                                     offset = offset,
                                     thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
          fitem$warning_lambda = lambda[iter]
          return(fitem)
        })
        if(!is.null(coef_res$warning_lambda)){
          warning_lambda = c(warning_lambda, coef_res$warning_lambda)
        }
        
        # Qui dobbiamo replicare il paradigma di active set
        
        # Non useremo aggiornare solo gli active ma tutti i beta dunque:
        #non Bmat[(Active[-c(1:k)] - k), iter] = coef_res$bvec[-(1:k)], ma:
        Bmat[ind_j[ind_j != (jth-k)],iter] = coef_res$bvec[-(1:k)]
        
        # Aggiungiamo una condizione di stop
        if (sum(Bmat[,iter] == 0) == length(Bmat[,iter])){
          break}
        
        # Di nuovo non consideriamo Active
        # Active <- as.numeric(c(1:k, sapply(which(Bmat[1:p, iter] != 0), function(i){
        #   k + seq(from = i, to = k * p, by = p)
        # })))
        
        # Nuovamente non aggiorniamo dunque i gruppi del prossimo fit:
        # group = c(1:k, rep((k+1):(k+length(which(Bmat[1:p, iter] != 0))), each = k))
        
        # Inoltre betaprev deve rimanere il risultato del fit
        # betaprev = c(coef_res$bvec[1:k], coef_res$bvec[-(1:k)][coef_res$bvec[-(1:k)]!=0]) # check intercette non nulle
        
        # anche qui annulliamo il betaprev
        # betaprev = coef_res$bvec
        
        
        b0[,iter] = coef_res$bvec[1:k]
        niterations[iter] <- coef_res$iterations  
      }
    }
  }
  return(list(b0 = b0, Bmat = Bmat, warning_lambda = warning_lambda, niterations = niterations))
}

# The fourth and last function is the inner EM algorithm that estimates NB parameters, given
# a variable "j", and given a lambda in the grid

zilgm_negbin2_grp = function(y, x, k, # dati e dimensioni
                             nk, inds, group, # indici e gruppi
                             theta = NULL, # dispersione
                             lambda, # penalty sparsità
                             bvec0 = NULL, # betahat del precedente lambda
                             w = NULL, penalty.factor = NULL, # penalty su n o p
                             offset = NULL,
                             # parametri tolleranza ciclo EM:
                             thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2, tol = 1e-6, maxit = 3e+2)
  {
  fun_call = match.call() # salva la chiamata alla funzioneper mostrarla in output
  
  # salviamo alcune variabili utili
  k = k
  N = NROW(x)
  P = NCOL(x) # 1500 quando privata delle jth-esime colonne
  p = (P - k )/k # era 499 quando chiamata con p = p-1
  
  if (!is.null(theta)) {
    fixed_theta = TRUE
    init_theta = theta
  } else {
    fixed_theta = FALSE
  }
  
  # Indici di elementi uguali a zero in y, necessario per fare inferenza su provenienza dalla delta o da NB
  pos_zero <- lapply(inds, function(i) y[i] == 0)
  pos_nzero <- lapply(pos_zero, function(x) !x)
  # Inizializzazione variabili latenti
  z <- lapply(nk, function(x) rep(1e-6, x))
  
  # Questa variabile la usiamo per la penalty lasssso, sarà la radice della cardinalità dei gruppi che moltiplica lambda, NB le intercette ne sono prive
  if (is.null(penalty.factor)) {
    penalty.factor <- c(rep(0, k), rep(sqrt(k), p))
  }
  
  if (is.null(w)) { # non è mai null in quanto inizializzato nella prima shell come tutti 1 comunque come check lo teniamo
    w = rep(1, N)
  }
  for(i in 1:k){ # deve essere 1/n_h per ogni h-esimo esperimento 
    w[inds[[i]]]= w[inds[[i]]]/length(inds[[i]]) 
  }# fattore di scala, si può elidere ed lasciare che sia "incorporato" su lambda se costante
  
  # Qui di seguito verifichiamo se i valori di y sono tutti uguali
  # In tal caso, cioè unique(y) == 1, il valore costante non potrebbe essere zero altrimenti avremmo filtrato quel gene in uno step di preprocessing
  # quindi qui si tiene solo conto del caso costante > 0 che comporterebbe prob = 0
  # come ulteriore conseguenza necessaria pos_zero sarà tutto falso
  
  if (length(unique(y)) == 1) {
    pos_zero <- lapply(pos_zero, function(x) rep(FALSE, length(x)))
    param = list(bvec = rep(0, P), sigma = 0, # bvec deve uscire 1500
                 # andrebbero inserite anche le intercette diverse da 0
                 prob = rep(0, k), pos_zero = pos_zero,
                 iterazio = 0)
    return(param)
  }
  
  # Adesso passiamo alla inizializzazione di parametri relativi al fit glm
  # Una stima iniziale del vettore delle medie:
  mu0 <- lapply(inds, function(i) {
    mu_ind <- y[i][y[i] > 0]
    rep(mean(mu_ind), length(i))
  })
  
  # Useremo beta in modo iterativo, innanzitutto potrebbe venir passato dalla stima del lambda precedente, qualora non fosse così però va inizializzato:
  if (is.null(bvec0)){
    bvec0 <- rep(0,(p+1)*k)
  } else {
    bvec0 = bvec0
  }
  # Passiamo adesso ad inizializzare il parametro che tiene conto della dispersione
  # Qui lo chiameremo theta come nella nostra formalizzazionema si riferisce in realtà all'alpha dell'articolo di park cioè il rapporto tra la media mu e il parametro di dispersione theta
  # con un theta alto questo rapporto deve risultare piccolo, la condizone iniziale è dunque simil poissoniana
  theta0 = 1e-4
  
  # Inizializziamo infine delle probabilità di avere zeri provenienti dalle delta
  prob0 <- sapply(1:k, function(i) {
    temp_prob <- (sum(pos_zero[[i]]) - sum(dNBII(0, mu = mu0[[i]], sigma = theta0, log = FALSE))) / nk[i]
    pmin(1, pmax(1e-10, temp_prob))
  })
  
  # Massimo della funzione obiettivo per la convergenza del ciclo EM
  erisk_prev = 1e+150
  
  
  iterazio = 0
  # Questo if apre al caso in cui nel vettore y non ci siano zeri, in tal caso non facciamo alcuna inferenza su prob, stimiamo direttamente media e dispersione
  if (sum(sapply(pos_zero, sum) == 0)) {
    ######################## CASO DA ADATTARE
    print("no zeroes in Y")
    # fit poisson per stimare beta
    sol_bvec = grpnetIACmu::grpnet.default(x = x, y = y, group = group, weights = w,
                                           family = "poisson", beta = bvec0,
                                           intercept = FALSE, lambda = lambda,
                                           standardize = FALSE, penalty = "LASSO",
                                           penalty.factor = penalty.factor,
                                           offset = offset,
                                           thresh = tol, maxit = maxit)
    
    nzeroel = sol_bvec$nzcoef # elementi non nulli
    bvec = sol_bvec$beta 
    mu = sol_bvec$mu
    # se il fit ha dato problemi di NaN nei coefficienti teniamo il beta precedente
    # Inoltre stimiamo theta con la formula inversa da media e varianza
    if (any(is.nan(bvec))){
      warning("NaN in the Poisson glm coefficients")
      bvec = bvec0
      theta = (mean(y)^2)/(var(y)-mean(y))
      prob = prob0
      erisk = erisk_prev
      nzeroel = sum(bvec0 != 0)
      mu = exp(offset + (x%*%bvec)/w)
    } else {
      # Altrimenti calcoliamo la nuova media dai beta
      # e poi la nuova dispersione
      #mu = exp(offset + (x%*%bvec)/w)
      if (fixed_theta) {
        theta = init_theta
      } else {
        theta = sigma_ml(y = y, mu = mu, weights = w)
      }
      # restituiamo l'output
      prob = prob0
      iterazio = 0
      erisk = erisk_prev
      theta = theta
      bvec0 = bvec  
    }
    
    ########################################
    # DOPO L'ELSE IL CICLO EM DI INTERESSE #
    ########################################
    
    
    
  } else {
    # Qui entriamo nel vero e proprio EM dove stimiamo iterativamente coefficienti, dispersione e probabilità di delta
    for (iterazio in 1:EM_iter) {
      print(iterazio)
      # Inizializzazione E-step: 
      prob = NULL
      tmp_z <- vector("list", k)
      
      # E-step: calcoliamo quanto sia probabile che uno zero venga dalla delta data l'attuale media mu0, il parametro theta0 (sempre rapporto tra mu0 e theta0 della NBI I vedi commento inizializzazione di theta0) e data l'attuale stima della prob
      # Il seguente codice sarebbe l'equazione 6 quando si usa pi invece che exp(fi)
      for (i in 1:k) {
        tmp_z[[i]] <- prob0[i] / (prob0[i] + (1 - prob0[i]) * dNBII(0, sigma = theta0, mu = mu0[[i]], log = FALSE))
        tmp_z[[i]][is.nan(tmp_z[[i]])] <- 1
        tmp_z[[i]] <- pmin(1 - 1e-6, tmp_z[[i]])
        
        z[[i]][pos_zero[[i]]] <- tmp_z[[i]][pos_zero[[i]]]
        
        prob[i] <- sum(z[[i]]) / nk[i]
        prob[i] <- pmin(1, pmax(1e-10, prob[i]))
      }
      ###############################################
      print(paste("le probabilità valgono:"))
      print(prob)
      ###############################################
      
      # M-step
      
      # Qui andiamo nello step di massimizzazione, partendo dal glm fit:
      # operiamo una grouped lasso poisson penalized regression
      # lo step precedente è tenuto in conto nei pesi w, che vengono aggiornati in modo tale da essere meno importanti nel fit quanto più è probabile che vengano dalla delta in zero
      
      sol_bvec = grpnetIACmu::grpnet.default(x = x, y = y, group = group, weights = w * (1 - unlist(z)),
                                             family = "poisson", beta = bvec0,
                                             intercept = FALSE, lambda = lambda, 
                                             standardize = FALSE, penalty = "LASSO",
                                             penalty.factor = penalty.factor,
                                             offset = offset,
                                             thresh = tol, maxit = maxit)
      
      # Estraiamo i risultati del fit
      
      nzeroel = sol_bvec$nzcoef # elementi non nulli
      bvec = sol_bvec$beta # coefficienti stimati (il supporto sono i nodi)
      mu = sol_bvec$mu
      # In particolare qui abbiamo 499*3 = 1497 coefficienti e 3 intercette
      # mancano i coefficienti di j vs j
      if (any(is.nan(bvec))){
        warning("NaN in the Poisson glm coefficients")
        bvec = bvec0
        theta = (mean(y)^2)/(var(y)-mean(y))
        prob = prob0
        erisk = erisk_prev
        nzeroel = sum(bvec0 != 0)
        mu = exp(offset + (x%*%bvec)/(w * (1 - unlist(z))))
        break
      }
      
      #mu = exp(offset + (x%*%bvec)/(w * (1 - unlist(z)))) # nuovo vettore delle medie
      
      # Qui se il coefficiente alpha era noto lo si tiene fissato altrimenti viene stimato partendo dalle medie risultanti dal fit e dall'attuale z
      
      if (fixed_theta) {
        theta = init_theta
      } else {
        theta = sigma_ml(y = y, mu = mu, weights = w * (1 - unlist(z)))
      }
      
      # Calcoliamo la funzione obiettivo
      erisk = grp_nb2_objective(y = y, prob = prob, bvec = bvec, mu = mu, lambda = lambda,
                                weights = w, penalty.factor = penalty.factor, sigma = theta,
                                posz = unlist(pos_zero), k = k, p = p)
      
      # Controlli di convergenza:
      # Se è infinito usa il valore precedente
      if (is.infinite(erisk) | is.nan(erisk)) {erisk = erisk_prev} 
      # nel caso in cui stiamo usando il valore precedente o non è cambiato entro una certa soglia fermiamo le iterazioni	  
      if ((abs((erisk_prev - erisk) / (erisk_prev + 1)) < EM_tol)) {
        bvec = bvec
        theta = theta
        prob = prob
        z = z
        break
        # } else if (erisk > erisk_prev + 1e-10) {
        #   bvec = bvec0
        #   theta = theta
        #   prob = prob0
        #   break
      } else { # se la funzione obiettivo può ancora crescere invece continua aggiornando i parametri:
        
        # erisk_prev = erisk
        # bvec0 = bvec
        # eta0 = eta
        # mu0 = mu
        # theta0 = theta
        # prob0 = prob
        
        erisk_prev = erisk
        mu0 = lapply(1:k, function(i){mu[inds[[i]]]})
        theta0 = theta
        prob0 = prob
        bvec0 = bvec
      }
    }
  }
  
  # flag = abs(bvec) < thresh
  # # stiamo annullando tutti i valori al di sotto di thresh!
  # bvec[flag] = 0
  # 
  # raccogliamo i risultati da inserire nell'output
  out = list()
  out$lambda = lambda
  out$bvec = bvec
  out$theta = theta
  out$prob = prob
  out$pos_zero = which(unlist(pos_zero))
  out$iterations = iterazio
  out$loglik = erisk
  out$call = fun_call
  out$nzeroel = nzeroel
  class(out) = "zilgm"
  return(out)
}
