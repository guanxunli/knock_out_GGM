library(doParallel)
library(RSpectra)
library(Seurat)
library(Matrix)
library(scTenifoldNet)

pcNet_ggm = function (X, nComp = 3, scaleScores = TRUE, symmetric = FALSE, 
                               q = 0, parallel = T, nclu = 4, verbose = TRUE) 
{
  if (!all(Matrix::rowSums(X) > 0)) {
    stop("Quality control has not been applied over the matrix.")
  }
  xClass <- class(X)[[1]]
  validClass <- xClass %in% c("matrix", "dgCMatrix")
  if (!validClass) {
    stop("Input should be a matrix with cells as columns and genes as rows")
  }
  if (nComp < 2 | nComp >= nrow(X)) {
    stop("nCom should be greater or equal than 2 and lower than the total number of genes")
  }
  gNames <- rownames(X)
  pcCoefficients <- function(K) {
    y <- X[, K]
    #Xi <- X
    Xi <- X[, -K]
    Xsvd = RSpectra::svds(Xi, nComp)
    coeff <- Xsvd$v
    Xsigma <- Xsvd$d
    score <- Xi %*% coeff
    score <- Matrix::t(Matrix::t(score)/(apply(score, 2, 
                                               function(X) {
                                                 sqrt(sum(X^2))
                                               })^2))
    Beta <- colSums(y * score) 
    ki_pre = 1 / var(c(y - score %*% Beta))
    Beta <- coeff %*% (Beta)
    return(c(Beta,ki_pre))
  }
  
  X <- (scale(Matrix::t(X)))
  n <- ncol(X)
  A <- 1 - diag(n)

    if (verbose) {
      B <- pbapply::pbsapply(seq_len(n), pcCoefficients)
    }
    else {
      B <- sapply(seq_len(n), pcCoefficients)
    }
  
  ki_pre = B[n,]
  B <- t(B[-n,])
  for (K in seq_len(n)) {
    A[K, A[K, ] == 1] = B[K, ]
  }
  
  if (isTRUE(symmetric)) {
    A <- (A + t(A))/2
  }
  absA <- abs(A)
  mA = max(absA)
  if (isTRUE(scaleScores)) {
    A <- (A/mA)
  }
  A[absA < quantile(absA, q)] <- 0
  diag(A) <- 1/ki_pre / mA
  colnames(A) <- rownames(A) <- gNames
  A <- as(A, "dgCMatrix")
  return(A)
}

makeNetworks_gmm = function (X, nNet = 10, nCells = 500, nComp = 3, scaleScores = TRUE, 
          symmetric = FALSE, q = 0.95, parallel = T, nclu = 4) 
{
  geneList <- rownames(X)
  nGenes <- length(geneList)
  nCol <- ncol(X)
  if (nGenes > 0) {
      pbapply::pbsapply(seq_len(nNet), function(W) {
        Z <- sample(x = seq_len(nCol), size = nCells, replace = TRUE)
        Z <- as.matrix(X[, Z])
        Z <- Z[apply(Z, 1, sum) > 0, ]
        if (nComp > 1 & nComp < nGenes) {
          Z <- pcNet_ggm(Z, nComp = nComp, scaleScores = scaleScores, 
                     symmetric = symmetric, q = q, verbose = FALSE)
        }
        else {
          stop("nComp should be greater or equal than 2 and lower than the total number of genes")
        }
        O <- matrix(data = 0, nrow = nGenes, ncol = nGenes)
        rownames(O) <- colnames(O) <- geneList
        O[rownames(Z), colnames(Z)] <- as.matrix(Z)
        O <- as(O, "dgCMatrix")
        return(O)
      },cl = nclu)
    }else {
    stop("Gene names are required")
  }
}


# add three new parameters
# 1) parallel: use parallel or not. If parallel = F, nclu is no use.
# 2) nclu: the number of nodes to use.
# 3) manifold_do: if manifold_do = F, do not do manifold alignment, only return two networks.
scTenifoldNet_gmm <- function(X, Y, qc_minLibSize = 1000, qc_removeOutlierCells = TRUE,
                           qc_minPCT = 0.05, qc_maxMTratio = 0.1, nc_nNet = 10,
                           nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                           nc_q = 0.05, td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5, parallel = F, nclu = NULL, ma_nDim = 30, dc_minFC = 1.5, manifold_do = F){
  # Single-cell Quality Control
  # X <- scQC(X, minLibSize = qc_minLibSize, removeOutlierCells = qc_removeOutlierCells, minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)
  # Y <- scQC(Y, minLibSize = qc_minLibSize, removeOutlierCells = qc_removeOutlierCells, minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)

  # Counts per million (CPM) normalization
  # X <- cpmNormalization(X)
  # Y <- cpmNormalization(Y)

  # Comparing gene ids.
  xNames <- rownames(X)
  yNames <- rownames(Y)
  
  sharedGenes <- intersect(xNames, yNames)
  nGenes <- length(sharedGenes)
  
  # Filtering out non-shared genes
  time1 = Sys.time()
  X <- X[sharedGenes,]
  Y <- Y[sharedGenes,]
  set.seed(1)
  xList <- makeNetworks_gmm(X = X, nCells = nc_nCells, nNet = nc_nNet, 
                        nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, 
                        q = (1 - nc_q), parallel = parallel, nclu = nclu)
  set.seed(2)
  yList <- makeNetworks_gmm(X = Y, nCells = nc_nCells, nNet = nc_nNet, 
                        nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, 
                        q = (1 - nc_q), parallel = parallel, nclu = nclu)
  xdiag = rep(0, nGenes)
  ydiag = rep(0, nGenes)
  for(i in 1:nc_nNet){
    xdiag = xdiag + 1/diag(xList[[i]])
    diag(xList[[i]]) = 0
    ydiag = ydiag + 1/diag(yList[[i]])
    diag(yList[[i]]) = 0
  }
  xdiag = xdiag / nc_nNet
  ydiag = ydiag / nc_nNet
  
  time2 = Sys.time()
  
  set.seed(1)
  tensorOut <- tensorDecomposition(xList, yList, K = td_K, 
                                   maxIter = td_maxIter, maxError = td_maxError)
  tX <- as.matrix(tensorOut$X)
  tY <- as.matrix(tensorOut$Y)
  
  diag(tX) = -1
  diag(tY) = -1
  tX = tX * (-xdiag)
  tY = tY * (-ydiag)
  
  set.seed(1)
  
  time3 = Sys.time()
  
  
  if(manifold_do){
    # 
    # Non-linear manifold alignment
    # for(A in c('O','D','P')){
    set.seed(1)
    mA <- manifoldAlignment(tX , tY, d = ma_nDim)
    rownames(mA) <- c(paste0('X_', sharedGenes),paste0('y_', sharedGenes))
    # outFile <-paste0(id,'_',M,'tensor_',A,'alignment.csv')
    # write.csv(mA, outFile)
    
    # Differential regulation testing
    dR <- dRegulation(manifoldOutput = mA, minFC = dc_minFC)
    # write.csv(dC, paste0('dCoex_',id,'_',M,'tensor_',A,'alignment.csv'))
    # }
    # }
    
    # Return preparation
    outputResult <- list()
    outputResult$Networks <- list(tX, tY)
    # outputResult$tensorNetworks <- tensorOut
    outputResult$manifoldAlignment <- mA
    outputResult$diffRegulation <- dR
  }else{
    outputResult <- list()
    outputResult$Networks <- list(tX, tY)
  }
  
  time4 = Sys.time()
  outputResult$time = c(time1, time2, time3, time4)
  
  # Return
  return(outputResult)
}


