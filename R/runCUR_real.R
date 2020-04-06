library(Matrix)
library(scTenifoldNet)
library(irlba)
library(rCUR) 

load('dataset/KO_real.RData')
WT <- as.matrix(real$tensorNetworks$X)

## read function and parameters
Args <- commandArgs()
k = as.numeric(Args[6])
n_rand = as.numeric(Args[7])
n_rand_c = as.numeric(Args[8])
n_rand_r = as.numeric(Args[9])
print(k)
print(n_rand_r)

## Do svd
print("Do svd.")
WT_svd <- irlba(WT, nv = k)


# source('R/utility.R')
gKO <- which(rownames(WT) %in% 'Nkx2-1')

## CUR original
CUR_ave_ori <- WT
CUR_ave_ori[gKO, ] <- 0
glist_cur_ori <- read.table("results/genelist/gList_cur_ori.txt")
glist_cur_ori <- glist_cur_ori$x
length(glist_cur_ori)

## CUR function
CUR_fun <- function(mat, mat_svd, gKO, k, n_rand, n_rand_c, n_rand_r){
  result_out <- list()
  result_sum <- 0
  for (i in 1:n_rand){
    CUR_res <- CUR(A = mat, c = n_rand_c, r = n_rand_r, sv = mat_svd, method = "exact.num.random")
    C <- CUR_res@C
    C[gKO, ] <- 0
    # C[, gKO] <- 0
    KT <- C %*% CUR_res@U %*% CUR_res@R
    result_out[[i]] <- KT
    result_sum <- result_sum + KT
  }
  result_sum <- result_sum / n_rand
  result_out[[n_rand + 1]] <-  result_sum
  return(result_out)
}

# run CUR function
print("run CUR function.")
CUR_res <- CUR_fun(mat = WT, mat_svd = WT_svd, gKO = gKO, k = k, n_rand = n_rand, 
                   n_rand_c = n_rand_c, n_rand_r = n_rand_r)

# save results
CUR_ave <- as.matrix(CUR_res[[n_rand + 1]])
colnames(CUR_ave) <- colnames(WT)
rownames(CUR_ave) <- rownames(WT)
saveRDS(CUR_ave, paste0("../results/CUR_ave_", k, "_", n_rand, "_", n_rand_c, ".rds"))

CUR_var <- matrix(0, nrow = nrow(CUR_ave), ncol = ncol(CUR_ave))
for (i in 1:n_rand){
  tmp <- as.matrix(CUR_res[[i]])
  CUR_var <- CUR_var + (tmp - CUR_ave)^2
}
CUR_var <- CUR_var / n_rand
saveRDS(CUR_var, paste0("../results/CUR_var_", k, "_", n_rand, "_", n_rand_c, ".rds"))

## run manifoldAlignment
print("Do manifold Alignment.")
cur_mA <- manifoldAlignment(WT, CUR_ave)
cur_DR <- dRegulation(cur_mA, minFC = 0)
save.image(paste0("CUR_real_", k, "_", n_rand_c, "_.Rdata"))

####################################
########## Analysis results ########
####################################

## 3_100_50
CUR_ave_3_100_50 <- readRDS("results/CUR_ave_3_100_50.rds")
sum((CUR_ave_3_100_50 - CUR_ave_ori)^2)

CUR_var_3_100_50 <- readRDS("results/CUR_var_3_100_50.rds")
range(CUR_var_3_100_50)

CUR_DR_3_100_50 <- readRDS("results/cur_DR_3_100_50.rds")
gList_cur_3_100_50 <- CUR_DR_3_100_50$gene[CUR_DR_3_100_50$p.value < 0.05]
length(gList_cur_3_100_50)
length(intersect(gList_cur_3_100_50, glist_cur_ori))

## 5_100_50
CUR_ave_5_100_50 <- readRDS("results/CUR_ave_5_100_50.rds")
sum((CUR_ave_5_100_50 - CUR_ave_ori)^2)
sum((CUR_ave_5_100_50 - CUR_ave_3_100_50)^2)

CUR_var_5_100_50 <- readRDS("results/CUR_var_5_100_50.rds")
range(CUR_var_5_100_50)

CUR_DR_5_100_50 <- readRDS("results/cur_DR_5_100_50.rds")
gList_cur_5_100_50 <- CUR_DR_5_100_50$gene[CUR_DR_5_100_50$p.value < 0.05]
length(intersect(gList_cur_5_100_50, glist_cur_ori))
length(intersect(gList_cur_5_100_50, gList_cur_3_100_50))
length(gList_cur_5_100_50)

## 10_100_50
CUR_ave_10_100_50 <- readRDS("results/CUR_ave_10_100_50.rds")
sum((CUR_ave_10_100_50 - CUR_ave_ori)^2)
sum((CUR_ave_10_100_50 - CUR_ave_3_100_50)^2)
sum((CUR_ave_10_100_50 - CUR_ave_5_100_50)^2)

CUR_var_10_100_50 <- readRDS("results/CUR_var_10_100_50.rds")
range(CUR_var_10_100_50)

CUR_DR_10_100_50 <- readRDS("results/cur_DR_10_100_50.rds")
gList_cur_10_100_50 <- CUR_DR_10_100_50$gene[CUR_DR_10_100_50$p.value < 0.05]
length(intersect(gList_cur_10_100_50, glist_cur_ori))
length(intersect(gList_cur_10_100_50, gList_cur_3_100_50))
length(intersect(gList_cur_10_100_50, gList_cur_5_100_50))
length(gList_cur_10_100_50)

## 3_100_100
CUR_ave_3_100_100 <- readRDS("results/CUR_ave_3_100_100.rds")
sum((CUR_ave_3_100_100 - CUR_ave_ori)^2)
sum((CUR_ave_3_100_100 - CUR_ave_3_100_50)^2)
sum((CUR_ave_3_100_100 - CUR_ave_5_100_50)^2)
sum((CUR_ave_3_100_100 - CUR_ave_10_100_50)^2)

CUR_var_3_100_100 <- readRDS("results/CUR_var_3_100_100.rds")
range(CUR_var_3_100_100)

CUR_DR_3_100_100 <- readRDS("results/cur_DR_3_100_100.rds")
gList_cur_3_100_100 <- CUR_DR_3_100_100$gene[CUR_DR_3_100_100$p.value < 0.05]
length(intersect(gList_cur_3_100_100, glist_cur_ori))
length(intersect(gList_cur_3_100_100, gList_cur_3_100_50))
length(intersect(gList_cur_3_100_100, gList_cur_5_100_50))
length(intersect(gList_cur_3_100_100, gList_cur_10_100_50))
length(gList_cur_3_100_100)

## 5_100_100
CUR_ave_5_100_100 <- readRDS("results/CUR_ave_5_100_100.rds")
sum((CUR_ave_5_100_100 - CUR_ave_ori)^2)
sum((CUR_ave_5_100_100 - CUR_ave_3_100_50)^2)
sum((CUR_ave_5_100_100 - CUR_ave_5_100_50)^2)
sum((CUR_ave_5_100_100 - CUR_ave_10_100_50)^2)
sum((CUR_ave_5_100_100 - CUR_ave_3_100_100)^2)

CUR_var_5_100_100 <- readRDS("results/CUR_var_5_100_100.rds")
range(CUR_var_5_100_100)

CUR_DR_5_100_100 <- readRDS("results/cur_DR_5_100_100.rds")
gList_cur_5_100_100 <- CUR_DR_5_100_100$gene[CUR_DR_5_100_100$p.value < 0.05]
length(intersect(gList_cur_5_100_100, glist_cur_ori))
length(intersect(gList_cur_5_100_100, gList_cur_3_100_50))
length(intersect(gList_cur_5_100_100, gList_cur_5_100_50))
length(intersect(gList_cur_5_100_100, gList_cur_10_100_50))
length(intersect(gList_cur_5_100_100, gList_cur_3_100_100))
length(gList_cur_5_100_100)

## 10_100_100
CUR_ave_10_100_100 <- readRDS("results/CUR_ave_10_100_100.rds")
sum((CUR_ave_10_100_100 - CUR_ave_ori)^2)
sum((CUR_ave_10_100_100 - CUR_ave_3_100_50)^2)
sum((CUR_ave_10_100_100 - CUR_ave_5_100_50)^2)
sum((CUR_ave_10_100_100 - CUR_ave_10_100_50)^2)
sum((CUR_ave_10_100_100 - CUR_ave_3_100_100)^2)
sum((CUR_ave_10_100_100 - CUR_ave_5_100_100)^2)

CUR_var_10_100_100 <- readRDS("results/CUR_var_10_100_100.rds")
range(CUR_var_10_100_100)

CUR_DR_10_100_100 <- readRDS("results/cur_DR_10_100_100.rds")
gList_cur_10_100_100 <- CUR_DR_10_100_100$gene[CUR_DR_10_100_100$p.value < 0.05]
length(intersect(gList_cur_10_100_100, glist_cur_ori))
length(intersect(gList_cur_10_100_100, gList_cur_3_100_50))
length(intersect(gList_cur_10_100_100, gList_cur_5_100_50))
length(intersect(gList_cur_10_100_100, gList_cur_10_100_50))
length(intersect(gList_cur_10_100_100, gList_cur_3_100_100))
length(intersect(gList_cur_10_100_100, gList_cur_5_100_100))
length(gList_cur_10_100_100)

