library(scTenifoldNet)
library(rCUR)
library(RSpectra)

load('dataset/KO_real.RData')
WT <- as.matrix(real$tensorNetworks$X)

gKO <- which(rownames(WT) %in% 'Nkx2-1')

load("results/GMM_CUR_real.Rdata")
length(gList_real)
length(gList_ggm_k)
length(intersect(gList_ggm_k, gList_real))

gList_link_row <- rownames(WT)[which(abs(WT[gKO, ]) > 0)]
length(gList_link_row)
length(intersect(gList_real, gList_link_row))
length(intersect(gList_ggm_k, gList_link_row))

tmp <- WT[gKO, which(abs(WT[gKO, ]) > 0)]
gList_link_row <- names(tmp)[order(abs(tmp), decreasing = TRUE)][1: (length(gList_real) * 3)]
length(intersect(gList_real, gList_link_row))
length(intersect(gList_ggm_k, gList_link_row))

gList_link_col <- rownames(WT)[which(abs(WT[, gKO]) > 0)]
length(gList_link_col)
length(intersect(gList_real, gList_link_col))
length(intersect(gList_ggm_k, gList_link_col))

tmp <- WT[which(abs(WT[, gKO]) > 0), gKO]
gList_link_col <- names(tmp)[order(abs(tmp), decreasing = TRUE)][1: (length(gList_real) * 3)]
length(intersect(gList_real, gList_link_col))
length(intersect(gList_ggm_k, gList_link_col))

dis_row <- rowSums((WT - matrix(WT[gKO, ], ncol = ncol(WT), nrow = nrow(WT), byrow = TRUE))^2)
gList_dis_row <- rownames(WT)[order(dis_row, decreasing = FALSE)][1: (length(gList_real) * 3)]
length(gList_dis_row)
length(intersect(gList_real, gList_dis_row))
length(intersect(gList_ggm_k, gList_dis_row))

dis_col <- colSums((WT - matrix(WT[, gKO], ncol = ncol(WT), nrow = nrow(WT)))^2)
gList_dis_col <- rownames(WT)[order(dis_col, decreasing = FALSE)][1: (length(gList_real) * 3)]
length(gList_dis_col)
length(intersect(gList_real, gList_dis_col))
length(intersect(gList_ggm_k, gList_dis_col))




