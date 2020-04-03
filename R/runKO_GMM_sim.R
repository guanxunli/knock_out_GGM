rm(list = ls())

library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(statsExpressions)
library(patchwork)
source("R/utility.R")

## input data
countMatrix <- read.csv('dataset/rn10g.csv', header = FALSE)
rownames(countMatrix) <- paste0('G', seq_len(nrow(countMatrix)))
countMatrix <- as.matrix(countMatrix)
countMatrix <- Matrix(countMatrix)

# logNormalization <- function(X){
#   index_rm <- which(colSums(X) == 0)
#   X <- X[, -index_rm]
#   X <- t(t(X)/colSums(X))*1e6
#   X <- as(X, 'dgCMatrix')
#   return(log(X + 1))
#   return(sqrt(X))
# }
# 
# countMatrix <- logNormalization(countMatrix)

## no difference
H0r <- runHO()

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

for(gKO in c(1,7)){
  print("Knock out at beginning, remove gene ", gKO)
  HAr <- runHA(gKO = gKO)
  D <- test_fun(H0r, HAr)
  print(D[order(D$A), ])
  
  # print("Knock out before tensor decomposition, remove gene ", gKO)
  # HAt <- runHAt(gKO = gKO)
  # D <- test_fun(H0r, HAt)
  # print(D[order(D$A), ])
  # 
  # print("Knock out before manifold alignment, remove gene ", gKO)
  # HAm <- runHAm(gKO = gKO)
  # D <- test_fun(H0r, HAm)
  # print(D[order(D$A), ])
  # 
  print("CUR decomposition original, remove gene ", gKO)
  HAc <- runHAc_v1(gKO = gKO)
  D <- test_fun(H0r, HAc)
  print(D[order(D$A), ])
  
  print("CUR decomposition random, remove gene ", gKO)
  HAc <- runHAc_v2(mat = countMatrix, gKO = gKO, n_rand = 100, n_rand_c = 5, n_rand_r = 5, k = 5)[[1]]
  D <- test_fun(H0r, HAc)
  print(D[order(D$A), ])
  
  print("GGM model w_ij = k_ij, remove gene ", gKO)
  Hggm1 <- runGGM_v1(gKO = gKO)
  D <- test_fun(H0r, Hggm1)
  print(D[order(D$A), ])
  
  print("GGM model w_ij = c_ij, remove gene ", gKO)
  Hggm2 <- runGGM_v2(gKO = gKO)
  D <- test_fun(H0r, Hggm2)
  print(D[order(D$A), ])
  
  print("GGM model w_ij = beta_ij, remove gene ", gKO)
  Hggm3 <- runGGM_v3(gKO = gKO)
  D <- test_fun(H0r, Hggm3)
  print(D[order(D$A), ])
}

source("R/scTenifold_subsample_ggm.R")
countMatrix <- as.matrix(countMatrix)
Hggm_y <- scTenifoldNet_gmm(countMatrix, countMatrix) #results are not good~ 