## Gaussian graphical model set w_ij as k_ij (should be symmetric)
GGM_v1 <- function(net_mat, index_rm, diag_net = 1){
  
  # set index
  n <- nrow(net_mat)
  index_per <- 1:n
  index_per <- index_per[-index_rm]
  
  # set diagnoal
  if (diag_net == 1){
    diag(net_mat) <- 1
  } else{
    diag(net_mat) <- sapply(1:n, FUN = function(x){return(max(c(net_mat[x, ], net_mat[, x])))}) * 1.1
  }
  
  ## get inverse of submatrix
  q11 <- net_mat[index_rm, index_rm]
  q1 <- net_mat[index_rm, -index_rm]
  q1t <- net_mat[-index_rm, index_rm]
  net_mat_inv <- net_mat[index_per, index_per] - tcrossprod(q1t, q1)/q11
  
  ## add matrix names
  if (is.null(rownames(net_mat))){
    return(net_mat_inv)
  } else{
    name_mat <- rownames(net_mat)
    rownames(net_mat_inv) <- name_mat[-index_rm]
    colnames(net_mat_inv) <- name_mat[-index_rm]
    return(net_mat_inv)
  }
}

## run GGM model with w_ij as k_ij
runGGM_v1 <- function(gKO, diag_net = "max"){
  WT <- countMatrix
  # Networks
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  
  # Tensor
  set.seed(1)
  WT <- tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  
  # GGM
  index <- which(rownames(WT) %in% paste0('G',gKO))
  temp <- GGM_v1(as.matrix(WT), index, diag_net = "max")
  
  KO <- matrix(0, nrow = nrow(WT), ncol = ncol(WT))
  rownames(KO) <- rownames(WT)
  colnames(KO) <- colnames(WT)
  diag(KO) <- diag(WT)
  KO[rownames(temp), colnames(temp)] <- temp
  KO[index, ] <- 0
  KO[, index] <- 0
  
  set.seed(1)
  mA <- manifoldAlignment(WT,KO, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}
