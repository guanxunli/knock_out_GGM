library(Rcpp)
library(SILGGM)
library(glasso)
library(igraph)
library(scTenifoldNet)
library(ggplot2)
library(GGally)
library(HurdleNormal)
library(dplyr)


## input data
countMatrix <- read.csv('python_SERGIO/simulationOutput.csv', header = FALSE)
rownames(countMatrix) <- paste0('G', seq_len(nrow(countMatrix)))
countMatrix <- as.matrix(countMatrix)

## network true
p <- ncol(countMatrix)
n <- nrow(countMatrix)
adj_mat <- matrix(0, ncol = n, nrow = n)
adj_mat[c(2,3,4,5,6), 1] <- 1
adj_mat[c(8,9,10), 7] <- 1
adj_mat_net <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
set.seed(2)
plot(adj_mat_net)

## Construct GGM
mat <- t(countMatrix)
res_B_NW_SL <- SILGGM(x = mat, method = "B_NW_SL")
res_D_S_NW_SL <- SILGGM(x = mat, method = "D-S_NW_SL")
res_D_S_GL <- SILGGM(x = mat, method = "D-S_GL")
res_GFC_SL <- SILGGM(x = mat, method = "GFC_SL")
# res_GFC_L <- SILGGM(x = mat, method = "GFC_L")

## check results

## B_NW_SL method
res_B_NW_SL$precision
res_B_NW_SL$partialCor
adj_precision <- matrix(0, ncol = n, nrow = n)
adj_precision[which(res_B_NW_SL$p_precision < 0.05)] <- 1
adj_precision_net <- graph_from_adjacency_matrix(adj_precision, mode = "undirected", diag = FALSE)
set.seed(2)
plot(adj_precision_net)

adj_partialCor <- matrix(0, ncol = n, nrow = n)
adj_partialCor[which(res_B_NW_SL$p_partialCor < 0.05)] <- 1
adj_partialCor_net <- graph_from_adjacency_matrix(adj_partialCor, mode = "undirected", diag = FALSE)
set.seed(2)
plot(adj_partialCor_net)

## D-S_NW_SL method failed
res_D_S_NW_SL$precision

## D-S_GL method
res_D_S_GL$precision
adj_precision <- matrix(0, ncol = n, nrow = n)
adj_precision[which(res_D_S_GL$p_precision < 0.05)] <- 1
adj_precision_net <- graph_from_adjacency_matrix(adj_precision, mode = "undirected", diag = FALSE)
set.seed(2)
plot(adj_precision_net)

## GFC_SL method
adj_test1 <- res_GFC_SL$global_decision[[1]]
adj_test_net <- graph_from_adjacency_matrix(adj_test1, mode = "undirected", diag = FALSE)
plot(adj_test_net)
adj_test2 <- res_GFC_SL$global_decision[[2]]
adj_test_net <- graph_from_adjacency_matrix(adj_test2, mode = "undirected", diag = FALSE)
set.seed(2)
plot(adj_test_net)


################################
######## My LASSO GGM ##########
################################
source("R/utility_GGM.R")
net_work <- GGM_net(countMatrix, method = "lasso", normlize = "FALSE", tk_qua = NULL, tk_rd = NULL)
adj_mat <- matrix(0, ncol = n, nrow = n)
adj_mat[which(abs(net_work) > 0)] <- 1
adj_mat_net <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
set.seed(2)
plot(adj_mat_net)

####################################################################
##############  ZERO-INFLATED SINGLE CELL GENE EXPRESSION ##########
####################################################################

dta_center <- conditionalCenter(mat) # 3000 x 10
hurdle_fit = fitHurdle(dta_center, fixed=NULL, parallel=FALSE, control=list(debug=0))
print(hurdle_fit)

## model selection based on BIC
BIC_etc <- hurdle_fit$BIC_etc
index_bic <- which.min(BIC_etc$BIC)
graph_seq <- BIC_etc[index_bic, ]
graphs_to_plot <- graph.adjacency(abs(hurdle_fit$adjMat[[graph_seq$adjMat_idx]])>0, mode='max')
hn_plot(graphs_to_plot, main=sprintf("Graph on %d edges", sum(degree(graphs_to_plot))/2))

## model selection based on true knowledge
ntureedges <- 8 * 2
index_test <- which.min((BIC_etc$trueEdges - ntureedges)^2)
graph_seq <- BIC_etc[index_test, ]
graphs_to_plot <- graph.adjacency(abs(hurdle_fit$adjMat[[graph_seq$adjMat_idx]])>0, mode='max')
hn_plot(graphs_to_plot, main=sprintf("Graph on %d edges", sum(degree(graphs_to_plot))/2))

################################
#### Recall scTenifoldnet ######
################################

## Adjacency matrix
WT <- countMatrix
set.seed(1)
WT <- makeNetworks(WT, nComp = 3, q = 0.95)
set.seed(1)
WT <- tensorDecomposition(WT)$X
WT <- as.matrix(WT)
WT <- (WT + t(WT))/2

adj_scTenifoldnet <- matrix(0, ncol = n, nrow = n)
adj_scTenifoldnet[which(abs(WT) > 0)] <- 1
adj_scTenifoldnet_net <- graph_from_adjacency_matrix(adj_scTenifoldnet, mode = "undirected", diag = FALSE)
set.seed(2)
plot(adj_scTenifoldnet_net)

## directed adjacency matrix
WT <- countMatrix
set.seed(1)
WT <- makeNetworks(WT, nComp = 3, q = 0.95)
set.seed(1)
WT <- tensorDecomposition(WT)$X
WT <- as.matrix(WT)

adj_scTenifoldnet_dir <- as.matrix(WT)
adj_scTenifoldnet_dir_net <- graph_from_adjacency_matrix(adj_scTenifoldnet_dir, mode = "directed", weighted = "TRUE", diag = FALSE)
diag(WT) <- 0
index <- which(abs(WT) > 0)
E(adj_scTenifoldnet_dir_net)$weight <- as.vector(abs(WT[index]))
E(adj_scTenifoldnet_dir_net)$width <- E(adj_scTenifoldnet_dir_net)$weight
edge_color <- ifelse(WT[index] > 0, "red", "blue")
E(adj_scTenifoldnet_dir_net)$edge.color <- edge_color
plot(adj_scTenifoldnet_dir_net, edge.arrow.size=.25, edge.color = edge_color)

