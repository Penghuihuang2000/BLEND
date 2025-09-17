BLEND_mixSQP <- function(Y, Ref, n.cores=NULL) {
  names.ref <- unlist(lapply(names(Ref),function(n){paste0(n,":",colnames(Ref[[n]]))}))
  ct.name <- names(Ref)
  Ref <- rlist::list.cbind(Ref)
  genes <- intersect(rownames(Y),rownames(Ref))
  Y <- Y[match(genes,rownames(Y)),]; Ref <- Ref[match(genes,rownames(Ref)),]
  
  cl <- parallel::makeCluster(ifelse(is.null(n.cores),max(parallel::detectCores()-1,1),n.cores))
  parallel::clusterExport(cl = cl, varlist = "Ref", envir = environment())
  out <- t(parallel::parApply(cl = cl, X = Y, MARGIN = 2, FUN = function(y){
    ind <- which(y>0)
    ashr::mixSQP(matrix_lik = Ref[ind,], weights = y[ind], prior = rep(1,ncol(Ref)))$pihat
  }))
  parallel::stopCluster(cl = cl)
  
  CTs <- sapply(strsplit(names.ref,":"),function(x){x[1]})
  subjs <- sapply(strsplit(names.ref,":"),function(x){x[2]})
  Pi.hat <- sapply(unique(CTs),function(ct){
    ind <- which(CTs==ct)
    if (length(ind)==1) {return(out[,ind])}
    return(rowSums(out[,ind]))
  }); colnames(Pi.hat) <- unique(CTs)
  if (!is.null(colnames(Y))) {rownames(Pi.hat) <- colnames(Y)}
  
  Mix.prop <- lapply(ct.name,function(ct){
    ind <- which(CTs==ct)
    if (length(ind)==1) {return(rep(1,NROW(Pi.hat)))}
    tmp.row.sum <- rowSums(out[,ind])
    tmp.row.sum <- sapply(tmp.row.sum, function(p){ifelse(p==0, 1e-15, p)})
    tmp <- out[,ind]/tmp.row.sum
    if (!is.null(colnames(Y))) {rownames(tmp) <- colnames(Y)}
    return(tmp)
  }); names(Mix.prop) <- unique(CTs)
  
  return(list(Pi.hat=Pi.hat,Mix.prop=Mix.prop))
}




