rm(list = ls())

library(igraph)
library(Matrix)
library(scTenifoldNet)

# Reading the input file
countMatrix <- read.csv('dataset/rn10g.csv', header = FALSE)
rownames(countMatrix) <- paste0('G', seq_len(nrow(countMatrix)))
countMatrix <- as.matrix(countMatrix)

# Input for scTenifoldNet
WT <- countMatrix

# Networks
set.seed(1)
WT <- makeNetworks(WT, nComp = 3, q = 0.8)
WT <- tensorDecomposition(WT)
WT <- as.matrix(WT$X)

source("R/utility.R")

## show the results
for(gKO in c(1,7)){
  print(paste0("Take absolute value, remove number ", gKO, " gene"))
  mkpositive_res <- mkpositive_fun(adj_mat = WT, KO_gene = gKO)
  print(mkpositive_res[order(abs(mkpositive_res), decreasing = TRUE)])
  
  print(paste0("Plus one, remove number ", gKO, " gene"))
  mkpositive_res <- mkpositive_fun(adj_mat = WT, KO_gene = gKO, method = "plusone")
  print(mkpositive_res[order(abs(mkpositive_res), decreasing = TRUE)])
  
  print(paste0("Remove weights, remove number ", gKO, " gene"))
  rmweight_res <- rmweight_fun(adj_mat = WT, KO_gene = gKO)
  print(rmweight_res[order(abs(rmweight_res), decreasing = TRUE)])
}


