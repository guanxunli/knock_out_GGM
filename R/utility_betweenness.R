library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(irlba)
library(igraph)

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
