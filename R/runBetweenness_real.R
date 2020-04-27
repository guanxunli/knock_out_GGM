library(scTenifoldNet)
library(rCUR)
library(patchwork)
library(igraph)

load('../dataset/KO_real.RData')
WT <- as.matrix(real$tensorNetworks$X)

source('R/utility_betweenness.R')
gKO <- which(rownames(WT) %in% 'Nkx2-1')

print("Do absolute value")
mkpositive_abs_res <- mkpositive_fun(adj_mat = WT, KO_gene = gKO)
print("Do plus one")
mkpositive_plus_res <- mkpositive_fun(adj_mat = WT, KO_gene = gKO, method = "plusone")
print("Do remove weights")
rmweight_res <- rmweight_fun(adj_mat = WT, KO_gene = gKO)

test_fun <- function(A){
  A <- as.data.frame(A)
  colnames(A) <- "Bet_infer"
  R <- sapply(seq_len(1e5), function(X){mean(sample(A$Bet_infer, replace = TRUE))})
  A$P <- sapply(A$Bet_infer, function(X){mean(R > X)})
  A$A <- p.adjust(A$P, method = 'fdr')
  A <- A[order(A$A), ]
  return(A)
}

mkpositive_abs_test <- test_fun(mkpositive_abs_res)
print(head(mkpositive_abs_test))
mkpositive_plus_test <- test_fun(mkpositive_plus_res)
print(head(mkpositive_plus_test))
rmweight_test <- test_fun(rmweight_res)
print(head(rmweight_test))


## check real results
source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
real$diffRegulation <- dRegulation(real$manifoldAlignment, minFC = 0)
gList_real <- real$diffRegulation$gene[real$diffRegulation$p.value < 0.05]
gList_abs <- rownames(mkpositive_abs_test)[mkpositive_abs_test$A < 0.05]
gList_plus <- rownames(mkpositive_plus_test)[mkpositive_plus_test$A < 0.05]
gList_rm <- rownames(rmweight_test)[rmweight_test$A < 0.05]

# save.image("Betweenness.Rdata")
load("results/Betweenness.Rdata")

print(length(gList_real))
print(length(gList_abs))
write.table(gList_abs, "gList_abs.txt")
print(length(intersect(gList_real, gList_abs)))
print(length(gList_plus))
write.table(gList_plus, "gList_plus.txt")
print(length(intersect(gList_real, gList_plus)))
print(length(gList_rm))
write.table(gList_rm, "gList_rm.txt")
print(length(intersect(gList_real, gList_rm)))

# mkpositive_abs_res <- mkpositive_abs_res[order(mkpositive_abs_res, decreasing = TRUE)]
# mkpositive_plus_res <- mkpositive_plus_res[order(mkpositive_plus_res, decreasing = TRUE)]
# rmweight_res <- rmweight_res[order(rmweight_res, decreasing = TRUE)]
# 
# length(intersect(gList_real, names(mkpositive_abs_res)[1:54]))
# length(intersect(gList_real, names(mkpositive_plus_res)[1:54]))
# length(intersect(gList_real, names(rmweight_res)[1:54]))



