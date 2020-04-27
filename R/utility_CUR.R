library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(statsExpressions)
library(patchwork)
library(irlba)
library(rCUR) 
library(RSpectra)

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
# gKO: The gene to be removed
# k: Truncated SVD
# n_rand: Number of tries
# n_rand_c: Number of column to be selected
# n_rand_r: Number of row to be selected

CUR_fun <- function(mat, gKO, k, n_rand, n_rand_c, n_rand_r){
  mat_svd <- irlba(mat, nv = k)
  result_out <- list()
  result_sum <- 0
  for (i in 1:n_rand){
    CUR_res <- CUR(A = mat, c = n_rand_c, r = n_rand_r, sv = mat_svd)
    C <- CUR_res@C
    C[gKO, ] <- 0
    # C[, gKO] <- 0
    KT <- C %*% CUR_res@U %*% CUR_res@R
    result_out[[i]] <- KT
    result_sum <- result_sum + KT
  }
  result_sum <- result_sum / n_rand
  result_out[[n_rand + 1]] <-  result_sum
  return(result_out)
}

runHAc_v2 <- function(mat, gKO, n_rand, n_rand_c, n_rand_r, k){
  WT <- mat
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  WT <- tensorDecomposition(WT)
  A <- WT$X
  A <- as.matrix(A)
  SVD <- RSpectra::svds(A, k = k)
  set.seed(1)
  curOutput <- CUR_fun(mat = A, gKO = gKO, k = k, n_rand = n_rand, n_rand_c = n_rand_c, n_rand_r = n_rand_r)
  B <- curOutput[[n_rand + 1]]
  colnames(B) <- colnames(A)
  rownames(B) <- rownames(A)
  set.seed(1)
  mA <- manifoldAlignment(A,B, d = 5)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste(1:10),]
  return(list("dR" = dR, "network" = curOutput))
}

#############################
########## test_fun #########
#############################

test_fun <- function(A, B){
  D <- reshape2::melt(list(WTvsWT = A$distance, WTvsKO = B$distance))
  D$Gene <- paste0('G', 1:10)
  D$Cluster <- c(rep('C1',6), rep('C2',4))
  colnames(D) <- c('Distance', 'Model', 'Gene', 'Cluster')
  D$Gene <- factor(D$Gene, levels = paste0('G', 1:10))
  D <- D[D$Model == 'WTvsKO',]
  D$FC <- B$FC
  D$Z <- scale(D$Distance)
  R <- sapply(seq_len(1e3), function(X){mean(sample(D$Z, replace = TRUE))})
  D$P <- sapply(D$Z, function(X){mean(R > X)})
  D$A <- p.adjust(D$P, method = 'fdr')
  D <- D[order(D$A, decreasing = TRUE), ]
  return(D)
}
