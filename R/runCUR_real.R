library(Matrix)
library(scTenifoldNet)
library(irlba)
library(rCUR) 

load('dataset/KO_real.RData')
WT <- as.matrix(real$tensorNetworks$X)
# WT_svd <- svds(WT, k = 100)

source('R/utility.R')
gKO <- which(rownames(WT) %in% 'Nkx2-1')

CUR_fun <- function(mat, gKO, k, n_rand, n_rand_c, n_rand_r){
  mat_svd <- irlba(mat, nv = k)
  result_out <- list()
  result_sum <- 0
  for (i in 1:n_rand){
    CUR_res <- CUR(A = mat, c = n_rand_c, r = n_rand_r, sv = mat_svd)
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

k <- 5
n_rand <- 100
n_rand_c <- 10
n_rand_d <- 10

CUR_res <- CUR_fun(mat = WT, gKO = gKO, k = k, n_rand = n_rand, n_rand_c = n_rand_c, n_rand_r = n_rand_d)
CUR_ave <- CUR_res[[n_rand + 1]]
colnames(CUR_ave) <- colnames(WT)
rownames(CUR_ave) <- rownames(WT)

cur_mA <- manifoldAlignment(WT, CUR_ave)
cur_DR <- dRegulation(cur_mA, minFC = 0)

save.image(paste0("CUR_real_", k, "_", n_rand_c, "_.Rdata"))