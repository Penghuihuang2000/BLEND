#' Title
#'
#' @import Rcpp
#' @import RcppArmadillo
#'
BLEND_EMMAP <- function(X_n, phi, alpha, beta, n.iter = 5000,
                     thres = 1e-5){
  # Dimensions
  G <- nrow(phi[[1]]) # number of genes
  CT <- length(phi) # number of cell types
  Ms <- unlist(lapply(phi, ncol)) # numbers of available references for cell types
  cumsum_Ms <- cumsum(Ms)
  cumsum_Ms_Ms <- cumsum_Ms - Ms + 1
  sourceCpp("src/code.cpp")
  # Initialization
  mu_n.i <- rep(1/CT, CT) # vector of length CT
  psi_n.i <- list() # list of length CT, each element has length Ms[[i]]
  for(i in 1:length(Ms)){
    psi_n.i <- c(psi_n.i, list(rep(1/Ms[i], Ms[i])))
  }
  U_n <- list() # list of length CT, each element is G X Ms[[i]]
  for(i in 1:length(Ms)){
    U_n <- c(U_n, list(matrix(0, nrow = G, ncol = Ms[i])))
  }
 check_list <- list()

  for(i in 1:n.iter){
    cat("\n iteration ",i)
    mu_n.old <- mu_n.i
    psi_n_mu_n <- list() # list of length CT, each element: vector Ms[ct]
    for(ct in 1:CT){
      psi_n_mu_n <- c(psi_n_mu_n, list(round(mu_n.i[ct]*psi_n.i[[ct]],6)))
    }

    result <- compute_U_n(X_n, phi, psi_n_mu_n, U_n,cumsum_Ms,cumsum_Ms_Ms,length(unlist(psi_n_mu_n)))
    U_n <- lapply(result, function(i) round(i,6))
    check_list[[i]] <- U_n
    for(ct in 1:CT){
      mu_n.i[ct] <- max(0, ((alpha+unlist(lapply(U_n, sum)))[ct]-1))
    }
    mu_n.i <- mu_n.i/sum(mu_n.i)
    for(ct in 1:CT){
      for(m in 1:(Ms[ct])){
        psi_n.i[[ct]][m] <- max(0, ((beta + colSums(U_n[[ct]]))[m] - 1))
      }
      if(sum(psi_n.i[[ct]])>0){
        psi_n.i[[ct]] <- psi_n.i[[ct]]/sum(psi_n.i[[ct]])
      }else{
        psi_n.i[[ct]] <- rep(1/length(psi_n.i[[ct]]), length(psi_n.i[ct]))
      }
    }
    if(sum(abs(mu_n.i - mu_n.old))<thres){
      return(list("cellular frac"=mu_n.i, "ref prop"=psi_n.i, "n.cvg"=i))
    }

    if((i %% 50) == 0) gc()
  }
 return(list("cellular frac"=mu_n.i, "ref prop"=psi_n.i, "n.cvg"=n.iter))
}


