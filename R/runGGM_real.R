library(scTenifoldNet)
library(rCUR)
library(patchwork)
library(RSpectra)

load('dataset/KO_real.RData')
WT <- as.matrix(real$tensorNetworks$X)
# WT_svd <- svds(WT, k = 100)

source('R/utility_GGM.R')
gKO <- which(rownames(WT) %in% 'Nkx2-1')

###############################################
########### GGM ###############################
##############################################

print("Do GGM")

# w_ij = k_ij
ggm_DR_k <- runGGM(net_mat = WT, gKO = gKO, diag_net = "max", wij = "kij", symme = "TRUE", rm0 = "FALSE")

# w_ij = k_ij
ggm_DR_c <- runGGM(net_mat = WT, gKO = gKO, diag_net = 1, wij = "cij", symme = "TRUE", rm0 = "FALSE")

# w_ij = k_ij
ggm_DR_b <- runGGM(net_mat = WT, gKO = gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "FALSE")

source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
real$diffRegulation <- dRegulation(real$manifoldAlignment, minFC = 0)

gList_real <- real$diffRegulation$gene[real$diffRegulation$p.value < 0.05]

gList_ggm_k <- ggm_DR_k$gene[ggm_DR_k$p.value < 0.05]
gList_ggm_c <- ggm_DR_c$gene[ggm_DR_c$p.value < 0.05]
gList_ggm_b <- ggm_DR_b$gene[ggm_DR_b$p.value < 0.05]

gList_cur_ori <- cur_DR_ori$gene[cur_DR_ori$p.value < 0.05]
gList_cur_rand <- cur_DR_rand$gene[cur_DR_rand$p.value < 0.05]


# save.image("GMM_CUR_real.Rdata")
load("results/GMM_CUR_real.Rdata")
length(gList_real)
write.table(gList_real, "gList_real.txt")
length(gList_ggm_k)
write.table(gList_ggm_k, "gList_ggm_k.txt")
length(intersect(gList_real, gList_ggm_k))
length(gList_ggm_c)
write.table(gList_ggm_c, "gList_ggm_c.txt")
length(intersect(gList_real, gList_ggm_c))
length(gList_ggm_b)
write.table(gList_ggm_b, "gList_ggm_b.txt")
length(intersect(gList_real, gList_ggm_b))
length(gList_cur_ori)
write.table(gList_cur_ori, "gList_cur_ori.txt")
length(intersect(gList_real, gList_cur_ori))
length(gList_cur_rand)
write.table(gList_cur_rand, "gList_cur_rand.txt")
length(intersect(gList_real, gList_cur_rand))


