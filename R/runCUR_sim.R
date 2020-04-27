rm(list = ls())

library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(statsExpressions)
library(patchwork)
source("R/utility_CUR.R")

## input data
countMatrix <- read.csv('dataset/rn10g.csv', header = FALSE)
rownames(countMatrix) <- paste0('G', seq_len(nrow(countMatrix)))
countMatrix <- as.matrix(countMatrix)
countMatrix <- Matrix(countMatrix)

## no difference
H0r <- runHO()

for(gKO in c(1,7)){
  print("Knock out at beginning, remove gene ", gKO)
  HAr <- runHA(gKO = gKO)
  D <- test_fun(H0r, HAr)
  print(D[order(D$A), ])
  
  print("Knock out before tensor decomposition, remove gene ", gKO)
  HAt<- runHAt(gKO = gKO)
  D <- test_fun(H0r, HAt)
  print(D[order(D$A), ])
  
  print("Knock out before manifold alignment, remove gene ", gKO)
  HAm <- runHAm(gKO = gKO)
  D <- test_fun(H0r, HAm)
  print(D[order(D$A), ])
  
  print("CUR original method, remove gene ", gKO)
  HAcur <- runHAc_v1(gKO = gKO)
  D <- test_fun(H0r, HAcur)
  print(D[order(D$A), ])
  
  print("CUR random method, remove gene ", gKO)
  HAcur <- runHAc_v2(gKO = gKO)
  D <- test_fun(H0r, HAcur)
  print(D[order(D$A), ])
}

