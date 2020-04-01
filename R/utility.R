library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(statsExpressions)
library(patchwork)

##############################################################
################## Basic try #################################
##############################################################

## Doesn't change anything
runHO <- function(){
  WT <- countMatrix
  KO <- countMatrix
  
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  KO <- makeNetworks(KO, nComp = 3, q = 0.8)
  
  set.seed(1)
  WT <- tensorDecomposition(WT)
  set.seed(1)
  KO <- tensorDecomposition(KO)
  
  set.seed(1)
  mA <- manifoldAlignment(WT$X,KO$X, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

## Knock out at the beginning
runHA <- function(gKO){
  WT <- countMatrix
  KO <- countMatrix
  KO[gKO,] <- 0
  
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  KO <- makeNetworks(KO, nComp = 3, q = 0.8)
  
  set.seed(1)
  WT <- tensorDecomposition(WT)
  set.seed(1)
  KO <- tensorDecomposition(KO)
  
  set.seed(1)
  mA <- manifoldAlignment(WT$X,KO$X, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

## knock out before tensor decomposition
runHAt <- function(gKO){
  WT <- countMatrix
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  KO <- WT
  for(i in 1:length(KO)){
    KO[[i]][gKO,] <- 0
    KO[[i]][,gKO] <- 0
  }
  
  set.seed(1)
  WT <- tensorDecomposition(WT)
  set.seed(1)
  KO <- tensorDecomposition(KO)
  
  set.seed(1)
  mA <- manifoldAlignment(WT$X,KO$X, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste(1:10),]
  return(dR)
}

## knock out before manifold alignment
runHAm <- function(gKO){
  WT <- countMatrix
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  WT <- tensorDecomposition(WT)
  KO <- WT
  KO$X[gKO,] <- 0
  KO$X[,gKO] <- 0
  set.seed(1)
  mA <- manifoldAlignment(WT$X,KO$X, d = 5)
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste(1:10),]
  return(dR)
}

###############################################
################## CUR ########################
###############################################

## CUR decomposition; Daniel's original method.
# runHAc_v1 <- function(gKO){
#   WT <- countMatrix
#   set.seed(1)
#   WT <- makeNetworks(WT, nComp = 3, q = 0.8)
#   set.seed(1)
#   WT <- tensorDecomposition(WT)
#   A <- WT$X
#   library(RSpectra)
#   SVD <- RSpectra::svds(A, k = 5)
#   library(rCUR)
#   set.seed(1)
#   curOutput <- rCUR::CUR(as.matrix(WT$X), sv = SVD)
#   curOutput@C[gKO,] <- 0 
#   #curOutput@C[,gKO] <- 0
#   #curOutput@U[gKO,] <- 0
#   #curOutput@U[,gKO] <- 0
#   #curOutput@R[gKO,] <- 0
#   #curOutput@R[,gKO] <- 0
#   B <- curOutput@C %*% curOutput@U %*% curOutput@R
#   colnames(B) <- colnames(A)
#   rownames(B) <- rownames(A)
#   set.seed(1)
#   mA <- manifoldAlignment(A,B, d = 5)
#   dR <- dRegulation(mA, minFC = 0)
#   dR <- dR[paste(1:10),]
#   return(dR)  
# }

## Daniel's version
runHAc_v1 <- function(gKO){
  WT <- countMatrix
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  WT <- tensorDecomposition(WT)
  A <- WT$X
  B <- A
  B[gKO, ] <- 0
  mA <- manifoldAlignment(A,B, d = 5)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste(1:10),]
  return(dR)  
}

## Try CUR
runHAc_v2 <- function(gKO, n_rand){
  WT <- countMatrix
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  WT <- tensorDecomposition(WT)
  A <- WT$X
  library(RSpectra)
  SVD <- RSpectra::svds(A, k = 5)
  library(rCUR)
  set.seed(1)
  curOutput <- rCUR::CUR(as.matrix(WT$X), sv = SVD)
  curOutput <- rCUR::CUR(as.matrix(WT$X), sv = SVD, c = n_rand, r = n_rand)
  curOutput@C[gKO,] <- 0 
  B <- curOutput@C %*% curOutput@U %*% curOutput@R
  colnames(B) <- colnames(A)
  rownames(B) <- rownames(A)
  set.seed(1)
  mA <- manifoldAlignment(A,B, d = 5)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste(1:10),]
  return(dR)  
}

##########################################################
################# GGM ####################################
##########################################################

## Gaussian graphical model set w_ij as k_ij (should be symmetric)
GGM_v1 <- function(net_mat, index_rm){
  
  n <- nrow(net_mat)
  if (index_rm == 1){
    per_net_mat <- net_mat
  } else if(index_rm == n){
    per_net_mat <- net_mat[c(index_rm, 1:(index_rm - 1)), c(index_rm, 1:(index_rm - 1))]
  } else{
    per_net_mat <- net_mat[c(index_rm, 1:(index_rm - 1), (index_rm + 1):n),
                           c(index_rm, 1:(index_rm - 1), (index_rm + 1):n)]
  }
  
  ## get inverse of submatrix
  q11 <- per_net_mat[1, 1]
  q1 <- per_net_mat[1, -1]
  q1t <- per_net_mat[-1, 1]
  net_mat_inv <- per_net_mat[2:n, 2:n] - tcrossprod(q1t, q1)/q11
  
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

## Gaussian graphical set w_ij as c_ij
GGM_v2 <- function(net_mat, index_rm){
  n <- nrow(net_mat)
  if (index_rm == 1){
    per_net_mat <- net_mat
  } else if(index_rm == n){
    per_net_mat <- net_mat[c(index_rm, 1:(index_rm - 1)), c(index_rm, 1:(index_rm - 1))]
  } else{
    per_net_mat <- net_mat[c(index_rm, 1:(index_rm - 1), (index_rm + 1):n),
                           c(index_rm, 1:(index_rm - 1), (index_rm + 1):n)]
  }
  
  ## get inverse of submatrix
  q11 <- per_net_mat[1, 1]
  q1 <- per_net_mat[1, -1]
  q1t <- per_net_mat[-1, 1]
  
  ## calculate factor
  tmp <- q1 * q1t
  fact_mat <- matrix(1, nrow = n - 1, ncol = n - 1) - matrix(tmp, nrow = n - 1, ncol = n - 1, byrow = TRUE) -
    matrix(tmp, nrow = n - 1, ncol = n - 1) + tcrossprod(tmp, tmp)
  net_mat_inv <- (per_net_mat[2:n, 2:n] + tcrossprod(q1t, q1)/q11)/sqrt(fact_mat)
  
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

## Gaussian graphical set w_ij as beta_ij
GGM_v3 <- function(net_mat, index_rm){
  n <- nrow(net_mat)
  if (index_rm == 1){
    per_net_mat <- net_mat
  } else if(index_rm == n){
    per_net_mat <- net_mat[c(index_rm, 1:(index_rm - 1)), c(index_rm, 1:(index_rm - 1))]
  } else{
    per_net_mat <- net_mat[c(index_rm, 1:(index_rm - 1), (index_rm + 1):n),
                           c(index_rm, 1:(index_rm - 1), (index_rm + 1):n)]
  }
  
  ## get inverse of submatrix
  # q11 <- per_net_mat[1, 1]
  q1 <- per_net_mat[1, -1]
  q1t <- per_net_mat[-1, 1]
  
  ## calculate factor
  tmp <- q1 * q1t
  fact_mat <- matrix(1, nrow = n - 1, ncol = n - 1) - matrix(tmp, nrow = n - 1, ncol = n - 1, byrow = TRUE)
  # net_mat_inv <- (per_net_mat[2:n, 2:n] + tcrossprod(q1t, q1)/q11)/fact_mat
  net_mat_inv <- (per_net_mat[2:n, 2:n] + tcrossprod(q1t, q1))/fact_mat
  
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
runGGM_v1 <- function(gKO){
  WT <- countMatrix
  # Networks
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  
  # Tensor
  set.seed(1)
  WT <- tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  diag(WT) <- 1
  # Your code
  index <- which(rownames(WT) %in% paste0('G',gKO))
  temp <- GGM_v1(as.matrix(WT), index)
  
  KO <- WT
  KO[KO != 0] <- 0
  diag(KO) <- 1
  KO[rownames(temp), colnames(temp)] <- temp
  
  set.seed(1)
  mA <- manifoldAlignment(WT,KO, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

## run GGM model with w_ij as c_ij
runGGM_v2 <- function(gKO){
  WT <- countMatrix
  # Networks
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  
  # Tensor
  set.seed(1)
  WT <- tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  diag(WT) <- 1
  # Your code
  WT <- (WT + t(WT))/2
  
  index <- which(rownames(WT) %in% paste0('G',gKO))
  temp <- GGM_v2(as.matrix(WT), index)
  
  KO <- WT
  KO[KO != 0] <- 0
  diag(KO) <- 1
  KO[rownames(temp), colnames(temp)] <- temp
  
  set.seed(1)
  mA <- manifoldAlignment(WT,KO, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

## run GGM model with w_ij as beta_ij
runGGM_v3 <- function(gKO){
  WT <- countMatrix
  # Networks
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  
  # Tensor
  set.seed(1)
  WT <- tensorDecomposition(WT)
  WT <- as.matrix(WT$X)
  # diag(WT) <- 1
  # Your code
  index <- which(rownames(WT) %in% paste0('G',gKO))
  temp <- GGM_v3(as.matrix(WT), index)
  
  KO <- WT
  KO[KO != 0] <- 0
  KO[rownames(temp), colnames(temp)] <- temp
  diag(KO) <- diag(WT)
  KO[gKO, gKO] <- 0
  KO <- round(KO, 1)
  
  
  set.seed(1)
  mA <- manifoldAlignment(WT,KO, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

###########################################################
#################### Betweenness ##########################
###########################################################

## make it positive weights
mkpositive_fun <- function(adj_mat, KO_gene, mode = "directed", weighted = TRUE, diag = FALSE, method = "abs"){
  
  if (method == "abs"){
    adj_mat <- abs(adj_mat)
  } else{
    adj_mat <- adj_mat + 1
  }
  
  graph_WT <- graph.adjacency(adj_mat, mode = mode, weighted = weighted, diag = diag)
  WT_betwenness <- betweenness(graph_WT, directed = (mode == "directed"))
  WT_r_betwenness <- WT_betwenness/sum(WT_betwenness)
  
  KO <- adj_mat[-KO_gene, -KO_gene]
  graph_KO <- graph.adjacency(KO, mode = mode, weighted = weighted, diag = diag)
  KO_betwenness <- betweenness(graph_KO, directed = (mode == "directed"))
  KO_r_betwenness <- KO_betwenness/sum(KO_betwenness)
  
  return(abs(WT_r_betwenness[-KO_gene] - KO_r_betwenness))
}

## remove weight infect
rmweight_fun <- function(adj_mat, KO_gene, threshold = 0.1, mode = "undirected", weighted = NULL, diag = FALSE){
  
  index0 <- (abs(adj_mat) < threshold)
  adj_mat_use <- matrix(1, nrow(adj_mat), ncol(adj_mat))
  adj_mat_use[index0] <- 0
  rownames(adj_mat_use) <- rownames(adj_mat)
  colnames(adj_mat_use) <- colnames(adj_mat)
  
  graph_WT <- graph.adjacency(adj_mat_use, mode = mode, weighted = weighted, diag = diag)
  WT_betwenness <- betweenness(graph_WT, directed = (mode == "directed"))
  WT_r_betwenness <- WT_betwenness/sum(WT_betwenness)
  
  KO <- adj_mat_use[-KO_gene, -KO_gene]
  graph_KO <- graph.adjacency(KO, mode = mode, weighted = weighted, diag = diag)
  KO_betwenness <- betweenness(graph_KO, directed = (mode == "directed"))
  KO_r_betwenness <- KO_betwenness/sum(KO_betwenness)
  
  return(abs(WT_r_betwenness[-KO_gene] - KO_r_betwenness))
}
