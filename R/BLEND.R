

#' @title BLEND
#' @description Estimate cell type fractions of bulk RNA-seq data using multiple references
#' @useDynLib BLEND
#' @param bulk bulk RNA-seq data. Counts. Genes by samples.
#' @param phi cell type-specific gene expression references. A list of length of the number of cell types. Each element of the list is a genes by reference matrix.
#' @param alpha numeric. prior for cellular fraction Dirichlet distribution.
#' @param beta numeric. prior for reference mixing proportion Dirichlet distribution.
#' @param ncore the number of cores that will be used for computation.
#' @param method posterior inference technique that will be used. Default is "mixSQP." Options are "mixSQP", "EMMAP", "GIBBS".
#' @param sample.idx vector, only for method = "GIBBS", index of Gibbs samples to be used for parameter estimation.
#' @param n.iter maximum iterations allowed in EM-MAP.
#' @param thres stopping criterion in EM-MAP.
#'
#' @return a list, cellular fractions and reference mixing proportions
#' @importFrom edgeR filterByExpr
#' @importFrom parallel makeCluster stopCluster parApply detectCores clusterExport
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom stats rgamma rmultinom
#' @importFrom rlist list.rbind list.cbind
#' @importFrom ashr mixSQP
#' @import Rcpp
#' @import RcppArmadillo
#' @export
#' 
BLEND <- function(bulk, phi, alpha = 1.00001, beta = 1.00001, ncore=50,
                  method = c("mixSQP", "EMMAP", "GIBBS"),
                  sample.idx=seq(from = 502, to = 2500, by = 2),
                  n.iter = 5000, thres = 1e-5){
  method <- match.arg(method)
  # Keep genes that are sufficient for statistical analysis
  # Very raw filtering
  keep.gene <- suppressMessages(edgeR::filterByExpr(bulk))
  keep.gene <- names(keep.gene)[keep.gene]
  bulk <- bulk[keep.gene,]
  
  # The situation where there's a cell type with only one reference
  ct_oneref <- c()
  for(i in 1:length(phi)){
    if(is.null(ncol(phi[[i]]))){
      ct_oneref <- c(ct_oneref, i)
      phi[[i]] <- as.matrix(cbind(phi[[i]], phi[[i]]))
      colnames(phi[[i]]) <- c("ref1","ref1_rep")
    }
  }

  # Record cell type names, bulk sample names, reference names
  # after computation, name outputs
  ct_name <- names(phi)
  if(is.null(colnames(bulk))){
    bulk_sample_name <- paste0("sample_", 1:ncol(bulk))
  }else{
    bulk_sample_name <- colnames(bulk)
  }
  reference_name.list <- list()
  for(i in 1:length(phi)){
    if(is.null(colnames(phi[[i]]))){
      reference_name <- paste0("reference_",1:ncol(phi[[i]]))
    }else{
      reference_name <- colnames(phi[[i]])
    }
    reference_name.list <- c(reference_name.list, list(reference_name))
  }
  names(reference_name.list) <- names(phi)
  
  
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
  if(! (method %in% c("mixSQP","GIBBS", "EMMAP"))){
    cat("\n Estimation approach not identified. \n Only mixSQP, EMMAP and GIBBS are supported!")
    return(0)
  }else{
    if(method == "GIBBS"){
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      est.res <- foreach(i = 1:ncol(bulk)) %dopar% {
        BLEND_GIBBS(X_n = bulk[,i], phi = phi, alpha = alpha, beta = beta,
                    sample.idx = sample.idx)
      }
      stopCluster(cl)
      est_frac <- rlist::list.rbind(lapply(est.res, function(x){x[[1]]}))
      est_ref_mix_prop <- list()
      for(i in 1:ncol(est_frac)){
        est_ref_mix_prop <- c(est_ref_mix_prop,
                              list(rlist::list.rbind(lapply(est.res, function(x){x[[2]][[i]]}))))
      }
    }
    if(method == "EMMAP"){
      cl <- makeCluster(ncore)
      registerDoParallel(cl)
      est.res <- foreach(i = 1:ncol(bulk), .errorhandling='pass', .packages = c("Rcpp", "RcppArmadillo","BLEND")) %dopar% {
        BLEND_EMMAP(X_n = bulk[,i], phi = phi, alpha = alpha, beta = beta,
                    n.iter = n.iter,
                    thres = thres)

      }
      stopCluster(cl)
      est_frac <- rlist::list.rbind(lapply(est.res, function(x){x[[1]]}))
      est_ref_mix_prop <- list()
      for(i in 1:ncol(est_frac)){
        est_ref_mix_prop <- c(est_ref_mix_prop,
                              list(rlist::list.rbind(lapply(est.res, function(x){x[[2]][[i]]}))))
      }
    }
    if(method == "mixSQP"){
      est.res <- BLEND_mixSQP(Y = bulk, Ref = phi, n.cores = ncore)
      est_frac <- est.res[[1]]
      est_ref_mix_prop <- est.res[[2]]
    }
  }
  gc()
  colnames(est_frac) <- ct_name
  rownames(est_frac) <- bulk_sample_name
  for(i in 1:length(est_ref_mix_prop)){
    rownames(est_ref_mix_prop[[i]]) <- bulk_sample_name
    colnames(est_ref_mix_prop[[i]]) <- reference_name.list[[i]]
  }
  names(est_ref_mix_prop) <- ct_name
  all_res <- list("cellular_fraction"=est_frac,
                  "reference_mix_prop"=est_ref_mix_prop)
  return(all_res)
}
