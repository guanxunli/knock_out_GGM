library(scTenifoldNet)
library(rCUR)
library(patchwork)
library(RSpectra)

load('dataset/KO_real.RData')
WT <- as.matrix(real$tensorNetworks$X)
WT_svd <- svds(WT, k = 100)

source('/R/utility.R')
gKO <- which(rownames(WT) %in% 'Nkx2-1')

##################################################
################### CUR ##########################
##################################################

print("Do CUR")
curKO_ori <- WT
curKO_ori[gKO,] <- 0

set.seed(1234)
curKO_rand <- CUR(WT, c = 100, r = 100, sv = WT_svd)
C <- curKO_rand@C
C[gKO,] <- 0
curKO_rand <- C %*% curKO_rand@U %*% curKO_rand@R

###############################################
########### GGM ###############################
##############################################

print("Do GGM")

# w_ij = k_ij
temp <- GGM_v1(WT, gKO)
ggmKO_k <- WT
ggmKO_k[ggmKO_k!=0] <- 0
ggmKO_k[rownames(temp), colnames(temp)] <- temp

# w_ij = c_ij
temp <- GGM_v2(WT, gKO)
ggmKO_c <- WT
ggmKO_c[ggmKO_c!=0] <- 0
ggmKO_c[rownames(temp), colnames(temp)] <- temp

# w_ij = beta_ij
temp <- GGM_v3(WT, gKO)
ggmKO_b <- WT
ggmKO_b[ggmKO_b!=0] <- 0
ggmKO_b[rownames(temp), colnames(temp)] <- temp


## manifold alignment
print("Do alignment")
set.seed(1)
ggm_mA_k <- manifoldAlignment(WT, ggmKO_k)
ggm_mA_c <- manifoldAlignment(WT, ggmKO_c)
ggm_mA_b <- manifoldAlignment(WT, ggmKO_b)
cur_mA_ori <- manifoldAlignment(WT, curKO_ori)
cur_mA_rand <- manifoldAlignment(WT, curKO_rand)

print("Do regulation")
ggm_DR_k <- dRegulation(ggm_mA_k, minFC = 0)
ggm_DR_c <- dRegulation(ggm_mA_c, minFC = 0)
ggm_DR_b <- dRegulation(ggm_mA_b, minFC = 0)
cur_DR_ori <- dRegulation(cur_mA_ori, minFC = 0)
cur_DR_rand <- dRegulation(cur_mA_rand, minFC = 0)

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


