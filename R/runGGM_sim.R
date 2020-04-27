rm(list = ls())

library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(statsExpressions)
library(patchwork)
source("R/utility_GGM.R")

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

##################################################################
################### construct network based on scTenifoldnet #####
##################################################################

for(gKO in c(1,7)){
  print("Knock out at beginning, remove gene ", gKO)
  HAr <- runHA(gKO = gKO)
  D <- test_fun(H0r, HAr)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symme = TRUE, rm0 = FALSE"))
  Hggm1 <- runGGM_sim(countMatrix, gKO, diag_net = "max", wij = "kij", symme = "TRUE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm1)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symcme = FALSE, rm0 = FALSE"))
  Hggm1 <- runGGM_sim(countMatrix, gKO, diag_net = "max", wij = "kij", symme = "c", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm1)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symme = TRUE, rm0 = TRUE"))
  Hggm1 <- runGGM_sim(countMatrix, gKO, diag_net = "max", wij = "kij", symme = "FALSE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm1, gKO = gKO)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = c_ij, remove gene ", gKO, " with symme = TRUE, rm0 = FALSE"))
  Hggm2 <- runGGM_sim(countMatrix, gKO = gKO, diag_net = 1, wij = "cij", symme = "TRUE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm2)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = c_ij, remove gene ", gKO, " with symme = TRUE, rm0 = TRUE"))
  Hggm2 <- runGGM_sim(countMatrix, gKO = gKO, diag_net = 1, wij = "cij", symme = "TRUE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm2, gKO = gKO)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = b_ij, remove gene ", gKO, " with symme = FALSE, rm0 = FALSE"))
  Hggm3 <- runGGM_sim(countMatrix, gKO = gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm3)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = b_ij, remove gene ", gKO, " with symme = FALSE, rm0 = TRUE"))
  Hggm3 <- runGGM_sim(countMatrix, gKO = gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm3, gKO = gKO)
  print(D[order(D$A), ])
}

###############################################################
########## construct network based on GGM #####################
###############################################################
# input data
countMatrix <- read.csv('dataset/rn10g.csv', header = FALSE)
rownames(countMatrix) <- paste0('G', seq_len(nrow(countMatrix)))
countMatrix <- as.matrix(countMatrix)
countMatrix <- Matrix(countMatrix)

net_work <- GGM_net(countMatrix, method = "gaussian", normlize = "FALSE", tk_qua = NULL, tk_rd = NULL)

colnames(net_work) <- paste0("G", 1:10)
rownames(net_work) <- paste0("G", 1:10)

## no difference
H0r <- runHO(WT = net_work)

gKO <- 1

for(gKO in c(1,7)){
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symme = TRUE, rm0 = FALSE"))
  Hggm1 <- runGGM(net_work, gKO, diag_net = "max", wij = "kij", symme = "TRUE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm1)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symcme = FALSE, rm0 = FALSE"))
  Hggm1 <- runGGM(net_work, gKO, diag_net = "max", wij = "kij", symme = "c", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm1)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = k_ij, remove gene ", gKO, " with symme = TRUE, rm0 = TRUE"))
  Hggm1 <- runGGM(net_work, gKO, diag_net = "max", wij = "kij", symme = "FALSE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm1, gKO = gKO)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = c_ij, remove gene ", gKO, " with symme = TRUE, rm0 = FALSE"))
  Hggm2 <- runGGM(net_work, gKO = gKO, diag_net = 1, wij = "cij", symme = "TRUE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm2)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = c_ij, remove gene ", gKO, " with symme = TRUE, rm0 = TRUE"))
  Hggm2 <- runGGM(net_work, gKO = gKO, diag_net = 1, wij = "cij", symme = "TRUE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm2, gKO = gKO)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = b_ij, remove gene ", gKO, " with symme = FALSE, rm0 = FALSE"))
  Hggm3 <- runGGM(net_work, gKO = gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "FALSE")
  D <- test_fun(H0r, Hggm3)
  print(D[order(D$A), ])
  
  print(paste0("GGM model w_ij = b_ij, remove gene ", gKO, " with symme = FALSE, rm0 = TRUE"))
  Hggm3 <- runGGM(net_work, gKO = gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "TRUE")
  D <- test_fun_rm0(H0r, Hggm3, gKO = gKO)
  print(D[order(D$A), ])
}

