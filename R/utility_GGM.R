library(Matrix)
library(scTenifoldNet)
library(ggplot2)
library(irlba)
library(rCUR) 
library(RSpectra)

##############################################################
################## Basic try #################################
##############################################################

## Doesn't change anything
runHO <- function(WT){
  set.seed(1)
  mA <- manifoldAlignment(WT,WT, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

##########################################################
################# GGM ####################################
##########################################################

## construct GGM net

GGM_net <- function(count_mat, method = "gaussian", normlize = "TRUE", tk_qua = NULL, tk_rd = NULL){
  n <- nrow(count_mat)
  net_work <- matrix(0, nrow = n, ncol = n)
  
  if (method == "gaussian"){
    for (i in 1:n){
      Y <- as.vector(countMatrix[i, ])
      X <- as.matrix(countMatrix[-i, ])
      lm.fit <- lm(Y ~ -1 + t(X))
      net_work[i, -i] <- coef(lm.fit)
    }
  } else if(method == "lasso"){
    for (i in 1:n){
      Y <- as.vector(countMatrix[i, ])
      X <- as.matrix(countMatrix[-i, ])
      glmnet.fit <- glmnet(x = t(X), y = Y, family = "gaussian", alpha = 1, intercept = FALSE,
                           lambda = cv.glmnet(x = t(X), y = Y)$lambda.1se)
      net_work[i, -i] <- coef(glmnet.fit)[-1, 1]
    }
  } else if(method == "glm"){
    for (i in 1:n){
      Y <- as.vector(countMatrix[i, ])
      X <- as.matrix(countMatrix[-i, ])
      glm.fit <- glm(Y ~ -1 + t(X), family = "poisson")
      net_work[i, -i] <- coef(glm.fit)
    }
  }
  
  # scale the network
  abs_net_work <- abs(net_work)
  if (normlize == TRUE){
    net_work <- net_work/max(abs_net_work)
  }
  
  # take round
  if(!is.null(tk_rd)){
    net_work <- round(net_work, tk_rd)
  }

  #remove small part
  if (!is.null(tk_qua)){
    net_work[abs_net_work < quantile(abs_net_work, tk_qua)] <- 0
  }
  
  return(net_work)
}

## GGM net after knock-out gene
GGM_invnet <- function(net_mat, index_rm, wij = "kij", symme = "TRUE"){
  
  n <- nrow(net_mat)
  index_per <- 1:n
  index_per <- c(index_rm, index_per[-index_rm])
  per_net_mat <- net_mat[index_per, index_per]
  
  ## get inverse of submatrix
  q11 <- per_net_mat[1, 1]
  q1 <- per_net_mat[1, -1]
  q1t <- per_net_mat[-1, 1]
  
  if(wij == "kij"){
    net_mat_inv <- per_net_mat[2:n, 2:n] - tcrossprod(q1t, q1)/q11
  } else if(wij == "cij"){
    tmp <- q1 * q1t
    fact_mat <- matrix(1, nrow = n - 1, ncol = n - 1) - matrix(tmp, nrow = n - 1, ncol = n - 1, byrow = TRUE) -
      matrix(tmp, nrow = n - 1, ncol = n - 1) + tcrossprod(tmp, tmp)
    net_mat_inv <- (per_net_mat[2:n, 2:n] - tcrossprod(q1t, q1))/sqrt(fact_mat)
  } else if(wij == "bij"){
    tmp <- q1 * q1t
    fact_mat <- matrix(1, nrow = n - 1, ncol = n - 1) - matrix(tmp, nrow = n - 1, ncol = n - 1)
    net_mat_inv <- (per_net_mat[2:n, 2:n] + tcrossprod(q1t, q1))/fact_mat
  }
  
  if (symme == "TRUE"){
    net_mat_inv <- (net_mat_inv + t(net_mat_inv))/2
  }
  
  ## add matrix names
  if (is.null(rownames(net_mat))){
    return(net_mat_inv)
  } else{
    name_mat <- rownames(net_mat)
    rownames(net_mat_inv) <- name_mat[-index_rm]
    colnames(net_mat_inv) <- name_mat[-index_rm]
    return(net_mat_inv)
  }
}


## run GGM model with w_ij as k_ij
runGGM <- function(net_mat, gKO, diag_net = "max", wij = "kij", symme = "TRUE", rm0 = "FALSE"){
  
  if (symme == "TRUE"){
    net_mat <- (net_mat + t(net_mat))/2
  }
  WT <- net_mat
  n <- nrow(WT)
  
  # set diagnoal
  if (diag_net == 1){
    diag(WT) <- 1
  } else{
    diag(WT) <- sapply(1:n, FUN = function(x){return(max(abs(c(WT[x, ], WT[, x]))))}) * 1.1
  }
  
  # index remove
  index_rm <- which(rownames(WT) %in% paste0('G',gKO))
  temp <- GGM_invnet(net_mat = as.matrix(WT), index_rm = index_rm, wij = wij, symme = symme)
  
  # getting new network
  KO <- matrix(0, nrow = nrow(WT), ncol = ncol(WT))
  rownames(KO) <- rownames(WT)
  colnames(KO) <- colnames(WT)
  
  # setting diagnoal
  if (wij == "bij" || wij == "cij"){
    diag(KO) <- diag(WT)
    # diag(KO) <- sapply(1:n, FUN = function(x){return(max(abs(c(KO[x, ], KO[, x]))))}) * 1.1
  }
  KO[rownames(temp), colnames(temp)] <- temp
  KO[index_rm, ] <- 0
  KO[, index_rm] <- 0
  
  if (rm0 =="TRUE"){
    KO <- KO[-index_rm, -index_rm]
  }
  
  # manifold alignment
  set.seed(1)
  mA <- manifoldAlignment(WT,KO, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

############################################
######## Modified GGm for simulation data ##
############################################

runGGM_sim <- function(net_mat, gKO, diag_net = "max", wij = "kij", symme = "TRUE", rm0 = "FALSE"){
  
  WT <- net_mat
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  WT <- tensorDecomposition(WT)$X
  WT <- as.matrix(WT)
  
  if (symme == "TRUE"){
    WT <- (WT + t(WT))/2
  }
  n <- nrow(WT)
  
  # set diagnoal
  if (diag_net == 1){
    diag(WT) <- 1
  } else{
    diag(WT) <- sapply(1:n, FUN = function(x){return(max(abs(c(WT[x, ], WT[, x]))))}) * 1.1
  }
  
  # index remove
  index_rm <- which(rownames(WT) %in% paste0('G',gKO))
  temp <- GGM_invnet(net_mat = as.matrix(WT), index_rm = index_rm, wij = wij, symme = symme)
  
  # getting new network
  KO <- matrix(0, nrow = nrow(WT), ncol = ncol(WT))
  rownames(KO) <- rownames(WT)
  colnames(KO) <- colnames(WT)
  
  # setting diagnoal
  if (wij == "bij" || wij == "cij"){
    diag(KO) <- diag(WT)
    # diag(KO) <- sapply(1:n, FUN = function(x){return(max(abs(c(KO[x, ], KO[, x]))))}) * 1.1
  }
  KO[rownames(temp), colnames(temp)] <- temp
  KO[index_rm, ] <- 0
  KO[, index_rm] <- 0
  
  if (rm0 =="TRUE"){
    KO <- KO[-index_rm, -index_rm]
  }
  
  # manifold alignment
  set.seed(1)
  mA <- manifoldAlignment(WT,KO, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

##########################################
######## test function ###################
##########################################
test_fun <- function(A, B){
  D <- reshape2::melt(list(WTvsWT = A$distance, WTvsKO = B$distance))
  D$Gene <- paste0('G', 1:10)
  D$Cluster <- c(rep('C1',6), rep('C2',4))
  colnames(D) <- c('Distance', 'Model', 'Gene', 'Cluster')
  D$Gene <- factor(D$Gene, levels = paste0('G', 1:10))
  D <- D[D$Model == 'WTvsKO',]
  D$FC <- B$FC
  D$Z <- scale(D$Distance)
  R <- sapply(seq_len(1e3), function(X){mean(sample(D$Z, replace = TRUE))})
  D$P <- sapply(D$Z, function(X){mean(R > X)})
  D$A <- p.adjust(D$P, method = 'fdr')
  D <- D[order(D$A, decreasing = TRUE), ]
  return(D)
}

test_fun_rm0 <- function(A, B, gKO){
  # set index
  n <- nrow(B)
  index_use <- 1:n
  index_use <- index_use[-gKO]
  A <- A[-gKO, ]
  B <- B[-n, ]
  # same as before
  D <- reshape2::melt(list(WTvsWT = A$distance, WTvsKO = B$distance))
  D$Gene <- paste0('G', index_use)
  colnames(D) <- c('Distance', 'Model', 'Gene')
  D$Gene <- factor(D$Gene, levels = paste0('G', index_use))
  D <- D[D$Model == 'WTvsKO',]
  D$FC <- B$FC
  D$Z <- scale(D$Distance)
  R <- sapply(seq_len(1e3), function(X){mean(sample(D$Z, replace = TRUE))})
  D$P <- sapply(D$Z, function(X){mean(R > X)})
  D$A <- p.adjust(D$P, method = 'fdr')
  D <- D[order(D$A, decreasing = TRUE), ]
  return(D)
}




