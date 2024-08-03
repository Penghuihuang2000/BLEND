## Generate sample from Dirichlet distribution given parameter
## Parameters:
### alpha: a real number vector
rdirichlet <- function(alpha){
  l <- length(alpha)
  x <- rgamma(length(alpha), alpha)
  if(identical(x, rep(0,l))){x[1] <- 1}
  return(x/sum(x))
}


## Estimate choice of reference, cellular fractions and updated
## individualized reference for a single bulk sample
## Parameters:
### X_n: counts vector of individual n, with length of G
### phi: a list of length CT, dimension G X M, Ms does not have to be the same
### alpha: prior for cellular fraction Dirichlet distribution
### beta: prior for reference proportion Dirichlet distribution
### sample.idx: index of Gibbs samples to be used for parameter estimation
BLEND_GIBBS <- function(X_n, phi, alpha, beta, sample.idx){
  # Dimensions
  G <- nrow(phi[[1]]) # number of genes
  CT <- length(phi) # number of cell types
  Ms <- unlist(lapply(phi, ncol)) # numbers of available references for cell types


  # Initialization
  mu_n.i <- rep(1/CT, CT) # vector of length CT
  psi_n.i <- list() # list of length CT, each element has length Ms[[i]]
  for(i in 1:length(Ms)){
    psi_n.i <- c(psi_n.i, list(rep(1/Ms[i], Ms[i])))
  }
  S_n <- list() # list of length CT, each element is G X Ms[[i]]
  for(i in 1:length(Ms)){
    S_n <- c(S_n, list(matrix(0, nrow = G, ncol = Ms[i])))
  }
  mu_n_sum <- rep(0, CT)
  psi_n_sum <- lapply(psi_n.i, function(x){0*x})
  S_n_sum <- S_n

  for(i in 1:max(sample.idx)){
    cat("\n iteration ",i)
    psi_n_mu_n <- list() # list of length CT, each element: vector Ms[ct]
    for(ct in 1:CT){
      psi_n_mu_n <- c(psi_n_mu_n, list(mu_n.i[ct]*psi_n.i[[ct]]))
    }
    for(g in 1:G){
      phi_g <- lapply(phi, function(x){x[g,]})
      prob.list <- psi_n_mu_n
      for(ct in 1:CT){
        prob.list[[ct]] <- prob.list[[ct]] * phi_g[[ct]]
      }
      S_n.sample <- rmultinom(n = 1,
                              size = X_n[g],
                              prob = unlist(prob.list))
      for(ct in 1:CT){
        S_n[[ct]][g,] <- S_n.sample[(cumsum(Ms) - Ms + 1)[ct]:(cumsum(Ms)[ct])]
      }
    }
    mu_n.i <- rdirichlet(alpha+unlist(lapply(S_n, sum)))
    for(ct in 1:CT){
      psi_n.i[[ct]] <- rdirichlet(beta + colSums(S_n[[ct]]))
    }

    if(i %in% sample.idx){
      mu_n_sum <- mu_n_sum + mu_n.i
      for(ct in 1:CT){
        psi_n_sum[[ct]] <- psi_n_sum[[ct]] + psi_n.i[[ct]]
        S_n_sum[[ct]] <- S_n_sum[[ct]] + S_n[[ct]]
      }
    }
    if((i %% 50) == 0) gc()
  }

  mu_n_est <- mu_n_sum/length(sample.idx)
  psi_n_est <- lapply(psi_n_sum, function(x){x/length(sample.idx)})
  return(list("cellular frac"=mu_n_est, "ref prop"=psi_n_est))

}


