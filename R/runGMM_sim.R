rm(list = ls())

library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(statsExpressions)
library(patchwork)
source("../R/utility.R")

## input data
countMatrix <- read.csv('../dataset/rn10g.csv', header = FALSE)
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

test_fun_rm0 <- function(A, B, gKO){
  # set index
  n <- nrow(B)
  index_use <- 1:n
  index_use <- index_use[-gKO]
  A <- A[-gKO, ]
  B <- B[-n, ]
  # same as before
  D <- reshape2::melt(list(WTvsWT = A$distance, WTvsKO = B$distance))
  D$Gene <- paste0('G', index_use)
  colnames(D) <- c('Distance', 'Model', 'Gene')
  D$Gene <- factor(D$Gene, levels = paste0('G', index_use))
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
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symme = TRUE, rm0 = FALSE"))
  Hggm1 <- runGGM(countMatrix, gKO, diag_net = "max", wij = "kij", symme = "TRUE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm1)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symcme = FALSE, rm0 = FALSE"))
  Hggm1 <- runGGM(countMatrix, gKO, diag_net = "max", wij = "kij", symme = "c", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm1)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symme = TRUE, rm0 = TRUE"))
  Hggm1 <- runGGM(countMatrix, gKO, diag_net = "max", wij = "kij", symme = "FALSE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm1, gKO = gKO)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = c_ij, remove gene ", gKO, " with symme = TRUE, rm0 = FALSE"))
  Hggm2 <- runGGM(countMatrix, gKO = gKO, diag_net = 1, wij = "cij", symme = "TRUE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm2)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = c_ij, remove gene ", gKO, " with symme = TRUE, rm0 = TRUE"))
  Hggm2 <- runGGM(countMatrix, gKO = gKO, diag_net = 1, wij = "cij", symme = "TRUE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm2, gKO = gKO)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = b_ij, remove gene ", gKO, " with symme = FALSE, rm0 = FALSE"))
  Hggm3 <- runGGM(countMatrix, gKO = gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm3)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = b_ij, remove gene ", gKO, " with symme = FALSE, rm0 = TRUE"))
  Hggm3 <- runGGM(countMatrix, gKO = gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm3, gKO = gKO)
  print(D[order(D$A), ])
}

