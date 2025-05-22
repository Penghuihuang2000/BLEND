

#' @title BLEND
#' @description Estimate cell type fractions of bulk RNA-seq data using multiple references
#' @useDynLib BLEND
#' @param bulk bulk RNA-seq data. Counts. Genes by samples.
#' @param phi cell type-specific gene expression references. A list of length of the number of cell types. Each element of the list is a genes by reference matrix.
#' @param alpha numeric. prior for cellular fraction Dirichlet distribution.
#' @param beta numeric. prior for reference mixing proportion Dirichlet distribution.
#' @param ncore the number of cores that will be used for computation.
#' @param method posterior inference technique that will be used. Default is "EMMAP." Can be chosen between "EMMAP" and "GIBBS".
#' @param sample.idx vector, index of Gibbs samples to be used for parameter estimation.
#' @param n.iter maximum iterations allowed in EM-MAP.
#' @param thres stopping criterion in EM-MAP.
#'
#' @return a list whose elements are the estimation results for each bulk sample
#' @importFrom edgeR filterByExpr
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom stats rgamma rmultinom
#' @export
#'
BLEND <- function(bulk, phi, alpha = 1.00001, beta = 1.00001, ncore=50,
                  method = c("EMMAP", "GIBBS"),
                  sample.idx=seq(from = 502, to = 2500, by = 2),
                  n.iter = 5000, thres = 1e-5){
  # Keep genes that are sufficient for statistical analysis
  # Very raw filtering
  keep.gene <- edgeR::filterByExpr(bulk)
  keep.gene <- names(keep.gene)[keep.gene]
  bulk <- bulk[keep.gene,]
  
  # The situation where there's a cell type with only one reference
  for(i in 1:length(phi)){
    if(is.null(ncol(phi[[i]]))){
      phi[[i]] <- as.matrix(cbind(phi[[i]], phi[[i]]))
      colnames(phi[[i]]) <- c("ref1","ref1_rep")
    }
  }

  # Bulk and reference shared genes
  reference.gene <- Reduce(intersect, lapply(phi, rownames))
  inter.gene <- intersect(rownames(bulk), reference.gene)
  if(length(inter.gene) == 0){
    cat("\n Genes no intersection! Please check gene names of bulk and reference.")
    return(0)
  }
  bulk <- bulk[inter.gene,]
  phi <- lapply(phi, function(x){x[inter.gene,]})

  # Drop genes that are not expressed in all cell types' references at all
  gene.names <- rownames(bulk)
  drop.gene <- c()
  for(i in 1:length(phi)){
    drop.gene <- c(drop.gene, which(rowSums(phi[[i]])==0))
  }
  drop.gene <- unique(drop.gene)
  if(length(drop.gene)!=0){
    bulk <- bulk[-drop.gene,]
    phi <- lapply(phi, function(x){x[-drop.gene,]})
    gene.names <- gene.names[-drop.gene]
  }

  # Normalize references
  phi <- lapply(phi, function(x){apply(x, 2, function(y){y/sum(y)})})

  # Use BLEND to deconvolve targeted bulk data
  library(foreach)
  library(doParallel)
  library(rlist)
  if(! (method %in% c("GIBBS", "EMMAP"))){
    cat("\n Estimation approach not identified. \n Only EMMAP and GIBBS are supported!")
    return(0)
  }else{
    if(method == "GIBBS"){
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      est.res <- foreach(i = 1:ncol(bulk)) %dopar% {
        # source("R/BLEND_GIBBS.R")
        BLEND_GIBBS(X_n = bulk[,i], phi = phi, alpha = alpha, beta = beta,
                    sample.idx = sample.idx)
      }
      stopCluster(cl)
    }
    if(method == "EMMAP"){
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      est.res <- foreach(i = 1:ncol(bulk), .errorhandling='pass', .packages = c("Rcpp", "RcppArmadillo","BLEND")) %dopar% {
        # source("./R/RcppExports.R")
        BLEND_EMMAP(X_n = bulk[,i], phi = phi, alpha = alpha, beta = beta,
                    n.iter = n.iter,
                    thres = thres)

      }
      stopCluster(cl)
    }
  }
  gc()
  return(est.res)
}
