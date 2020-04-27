## load data
rm(list = ls())

library(Matrix)
library(scTenifoldNet)
library(ggplot2)

## manifold alignment
manifoldAlignment <- function(X, Y, d = 30){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X+1
  wY <- Y+1
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- E$vectors[,seq_len(d)]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}

## Doesn't change anything
runHO <- function(){
  WT <- countMatrix
  KO <- countMatrix
  
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  KO <- makeNetworks(KO, nComp = 3, q = 0.8)
  
  set.seed(1)
  WT <- tensorDecomposition(WT)
  set.seed(1)
  KO <- tensorDecomposition(KO)
  
  set.seed(1)
  mA <- manifoldAlignment(WT$X,KO$X, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(dR)
}

## Daniel's version
runHAc_v1 <- function(gKO){
  WT <- countMatrix
  set.seed(1)
  WT <- makeNetworks(WT, nComp = 3, q = 0.8)
  set.seed(1)
  WT <- tensorDecomposition(WT)
  A <- WT$X
  B <- A
  B[gKO, ] <- 0
  mA <- scTenifoldNet::manifoldAlignment(A,B, d = 5)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste(1:10),]
  return(dR)  
}


## test function
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

## make network
makeNetworks <- function (X, nNet = 10, nCells = 500, nComp = 3, scaleScores = TRUE, 
                          symmetric = FALSE, q = 0.95) 
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
        Z <- pcNet(Z, nComp = nComp, scaleScores = scaleScores, 
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
    })
  }
  else {
    stop("Gene names are required")
  }
}

pcNet <- function(X, nComp = 3, scaleScores = TRUE, symmetric = FALSE, 
                   q = 0, verbose = TRUE) {
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
    Xi <- X
    Xi <- Xi[, -K]
    coeff <- RSpectra::svds(Xi, nComp)$v
    score <- Xi %*% coeff
    score <- Matrix::t(Matrix::t(score)/(apply(score, 2, 
                                               function(X) {
                                                 sqrt(sum(X^2))
                                               })^2))
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)
    return(Beta)
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
  B <- t(B) # every row is one result for one regression
  for (K in seq_len(n)) {
    A[K, A[K, ] == 1] = B[K, ]
  }
  if (isTRUE(symmetric)) {
    A <- (A + t(A))/2
  }
  absA <- abs(A)
  if (isTRUE(scaleScores)) {
    A <- (A/max(absA))
  }
  A[absA < quantile(absA, q)] <- 0
  diag(A) <- 0
  colnames(A) <- rownames(A) <- gNames
  A <- as(A, "dgCMatrix")
  return(A)
}

tensorDecomposition <- function (xList, yList = NULL, K = 5, maxError = 1e-05, maxIter = 1000) 
{
  xNets <- length(xList)
  if (!is.null(yList)) {
    yNets <- length(yList)
    if (xNets != yNets) {
      stop("Same number of networks are required in both cases")
    }
    nNet <- unique(c(xNets, yNets))
    xGenes <- unique(unlist(lapply(xList, rownames)))
    yGenes <- unique(unlist(lapply(yList, rownames)))
    sGenes <- intersect(xGenes, yGenes)
  }
  else {
    nNet <- xNets
    xGenes <- unique(unlist(lapply(xList, rownames)))
    sGenes <- xGenes
  }
  nGenes <- length(sGenes)
  tensorX <- array(data = 0, dim = c(nGenes, nGenes, 1, nNet))
  if (!is.null(yList)) {
    tensorY <- array(data = 0, dim = c(nGenes, nGenes, 1, 
                                       nNet))
  }
  for (i in seq_len(nNet)) {
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes, tGenes] <- temp[tGenes, tGenes]
    tensorX[, , , i] <- tempX
    if (!is.null(yList)) {
      tempY <- matrix(0, nGenes, nGenes)
      rownames(tempY) <- colnames(tempY) <- sGenes
      temp <- as.matrix(yList[[i]])
      tGenes <- sGenes[sGenes %in% rownames(temp)]
      tempY[tGenes, tGenes] <- temp[tGenes, tGenes]
      tensorY[, , , i] <- tempY
    }
  }
  set.seed(1)
  tensorX <- as.tensor(tensorX)
  tensorX <- cpDecomposition(tnsr = tensorX, num_components = K, 
                             max_iter = maxIter, tol = maxError)
  tX <- tensorX$est$data[, , , 1]
  for (i in seq_len(nNet)[-1]) {
    tX <- tX + tensorX$est$data[, , , i]
  }
  tX <- tX/nNet
  # tX <- tX/max(abs(tX))
  # tX <- round(tX, 1)
  tX <- as(tX, "dgCMatrix")
  rownames(tX) <- colnames(tX) <- sGenes
  if (!is.null(yList)) {
    set.seed(1)
    tensorY <- as.tensor(tensorY)
    tensorY <- cpDecomposition(tnsr = tensorY, num_components = K, 
                               max_iter = 1000)
    tY <- tensorY$est$data[, , , 1]
    for (i in seq_len(nNet)[-1]) {
      tY <- tY + tensorY$est$data[, , , i]
    }
    tY <- tY/nNet
    # tY <- tY/max(abs(tY))
    # tY <- round(tY, 1)
    tY <- as(tY, "dgCMatrix")
    rownames(tY) <- colnames(tY) <- sGenes
  }
  tensorOutput <- list()
  tensorOutput$X <- tX
  if (!is.null(yList)) {
    tensorOutput$Y <- tY
  }
  return(tensorOutput)
}

cpDecomposition <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5){
  kronecker_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for(i in 2:length(L)){
      retmat <- kronecker(retmat,L[[i]])
    }
    retmat
  }
  
  
  fnorm <- function(tnsr){
    arr<-tnsr$data
    sqrt(sum(arr*arr))
  }
  
  rs_unfold <- function(tnsr,m=NULL){
    if(is.null(m)) stop("mode m must be specified")
    num_modes <- tnsr$num_modes
    rs <- m
    cs <- (1:num_modes)[-m]
    unfold(tnsr,row_idx=rs,col_idx=cs)
  }
  
  unfold <- function(tnsr,row_idx=NULL,col_idx=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row and column indices must be specified")
    num_modes <- tnsr$num_modes
    if (length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
    if(any(rs<1) || any(rs>num_modes) || any(cs < 1) || any(cs>num_modes)) stop("illegal indices specified")
    perm <- c(rs,cs)
    if (any(sort(perm,decreasing=TRUE) != num_modes:1)) stop("missing and/or repeated indices")
    modes <- tnsr$modes
    mat <- tnsr$data
    new_modes <- c(prod(modes[rs]),prod(modes[cs]))
    #rearranges into a matrix
    mat <- aperm(mat,perm)
    dim(mat) <- new_modes
    as.tensor(mat)
  }
  
  hadamard_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for (i in 2:length(L)){
      retmat <- retmat*L[[i]]
    }
    retmat
  }
  
  khatri_rao_list <- function(L,reverse=FALSE){
    stopifnot(all(unlist(lapply(L,is.matrix))))
    ncols <- unlist(lapply(L,ncol))
    stopifnot(length(unique(ncols))==1)
    ncols <- ncols[1]
    nrows <- unlist(lapply(L,nrow))
    retmat <- matrix(0,nrow=prod(nrows),ncol=ncols)
    if (reverse) L <- rev(L)
    for(j in 1:ncols){
      Lj <- lapply(L,function(x) x[,j])
      retmat[,j] <- kronecker_list(Lj)
    }
    retmat
  }
  
  superdiagonal_tensor <- function(num_modes,len,elements=1L){
    modes <- rep(len,num_modes)
    arr <- array(0, dim = modes)
    if(length(elements)==1) elements <- rep(elements,len)
    for (i in 1:len){
      txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
      eval(parse(text=txt))
    }
    as.tensor(arr)
  }
  
  ttl<-function(tnsr,list_mat,ms=NULL){
    if(is.null(ms)||!is.vector(ms)) stop ("m modes must be specified as a vector")
    if(length(ms)!=length(list_mat)) stop("m modes length does not match list_mat length")
    num_mats <- length(list_mat)
    if(length(unique(ms))!=num_mats) warning("consider pre-multiplying matrices for the same m for speed")
    mat_nrows <- vector("list", num_mats)
    mat_ncols <- vector("list", num_mats)
    for(i in 1:num_mats){
      mat <- list_mat[[i]]
      m <- ms[i]
      mat_dims <- dim(mat)
      modes_in <- tnsr$modes
      stopifnot(modes_in[m]==mat_dims[2])
      modes_out <- modes_in
      modes_out[m] <- mat_dims[1]
      tnsr_m <- rs_unfold(tnsr,m=m)$data
      retarr_m <- mat%*%tnsr_m
      tnsr <- rs_fold(retarr_m,m=m,modes=modes_out)
    }	
    tnsr	
  }
  
  rs_fold <- function(mat,m=NULL,modes=NULL){
    if(is.null(m)) stop("mode m must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    num_modes <- length(modes)
    rs <- m
    cs <- (1:num_modes)[-m]
    fold(mat,row_idx=rs,col_idx=cs,modes=modes)
  }
  
  fold <- function(mat, row_idx = NULL, col_idx = NULL, modes=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    if(!is(mat,"list")){
      if(!is.matrix(mat))  stop("mat must be of class 'matrix'")
    }else{
      stopifnot(mat$num_modes==2)
      mat <- mat$data			
    }
    num_modes <- length(modes)
    stopifnot(num_modes==length(rs)+length(cs))
    mat_modes <- dim(mat)
    if((mat_modes[1]!=prod(modes[rs])) || (mat_modes[2]!=prod(modes[cs]))) stop("matrix nrow/ncol does not match Tensor modes")
    #rearranges into array
    iperm <- match(1:num_modes,c(rs,cs))
    as.tensor(aperm(array(mat,dim=c(modes[rs],modes[cs])),iperm))
  }
  
  
  if(is.null(num_components)) stop("num_components must be specified")
  stopifnot(is(tnsr,"list"))
  #if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
  
  #initialization via truncated hosvd
  num_modes <- tnsr$num_modes
  modes <- tnsr$modes
  U_list <- vector("list",num_modes)
  unfolded_mat <- vector("list",num_modes)
  tnsr_norm <- fnorm(tnsr)
  for(m in 1:num_modes){
    unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)$data
    U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
  }
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  #set up convergence check
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(est){
    curr_resid <- fnorm(as.tensor(est$data - tnsr$data))
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter==1) return(FALSE)
    if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
    else{ return(FALSE)}
  }	
  #progress bar
  pb <- txtProgressBar(min=0,max=max_iter,style=3)
  #main loop (until convergence or max_iter)
  norm_vec <- function(vec){
    Matrix::norm(as.matrix(vec))
  }
  while((curr_iter < max_iter) && (!converged)){
    setTxtProgressBar(pb,curr_iter)
    for(m in 1:num_modes){
      V <- hadamard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
      V_inv <- solve(V)			
      tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
      lambdas <- apply(tmp,2,norm_vec)
      U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
      Z <- superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
      est <- ttl(Z,U_list,ms=1:num_modes)
    }
    #checks convergence
    if(CHECK_CONV(est)){
      converged <- TRUE
      setTxtProgressBar(pb,max_iter)
    }else{
      curr_iter <- curr_iter + 1
    }
  }
  if(!converged){setTxtProgressBar(pb,max_iter)}
  close(pb)
  #end of main loop
  #put together return list, and returns
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]
  norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
  invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
}

as.tensor <- function(x,drop=FALSE){
  stopifnot(is.array(x)||is.vector(x))
  tnsr <- list()
  if (is.vector(x)){
    modes <- c(length(x))
    num_modes <- 1L
    tnsr$modes <- modes
    tnsr$num_modes <- num_modes
    tnsr$data <- x
    #new("Tensor", num_modes, modes, data = x)
  }
  else {
    modes <- dim(x)
    num_modes <- length(modes)
    dim1s <- which(modes==1)
    if (drop && (length(dim1s)>0)){
      modes <- modes[-dim1s]
      num_modes <- num_modes-length(dim1s)
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- array(x,dim=modes)
      #new("list",num_modes,modes,data=array(x,dim=modes))
    } else {
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- x
    }
  }
  return(tnsr)
}

## GGM net
GGM_invnet <- function(net_mat, index_rm, wij = "kij", symme = "TRUE"){
  
  net_mat <- as.matrix(net_mat)
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
  
  net_mat <- as.matrix(net_mat)
  if (symme == "TRUE"){
    net_mat <- (net_mat + t(net_mat))/2
  }
  WT <- net_mat
  n <- nrow(net_mat)
  
  # set diagnoal
  if (diag_net == 1){
    diag(net_mat) <- 1
  } else{
    diag(net_mat) <- sapply(1:n, FUN = function(x){return(max(abs(c(WT[x, ], WT[, x]))))}) * 1.1
  }
  
  # index remove
  index_rm <- which(rownames(net_mat) %in% paste0('G',gKO))
  temp <- GGM_invnet(net_mat = as.matrix(net_mat), index_rm = index_rm, wij = wij, symme = symme)
  
  # getting new network
  KO <- matrix(0, nrow = nrow(net_mat), ncol = ncol(net_mat))
  rownames(KO) <- rownames(net_mat)
  colnames(KO) <- colnames(net_mat)
  
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
  
  diag(WT) <- 0
  diag(KO) <- 0
  
  mA <- manifoldAlignment(WT,KO, d = 5)
  
  set.seed(1)
  dR <- dRegulation(mA, minFC = 0)
  dR <- dR[paste0(1:10),]
  return(list(dR, WT, KO))
}


## input data
countMatrix <- read.csv('dataset/rn10g.csv', header = FALSE)
rownames(countMatrix) <- paste0('G', seq_len(nrow(countMatrix)))
countMatrix <- as.matrix(countMatrix)
countMatrix <- Matrix(countMatrix)

WT <- countMatrix
WT <- makeNetworks(WT, nNet = 10, nCells = 1000, nComp = 3, scaleScores = FALSE, 
                   symmetric = FALSE, q = 0.8)
WT <- tensorDecomposition(WT)$X

## no difference
H0r <- runHO()
H0r_new <- runHO()
## original method
gKO <- 1
Hggm <- runGGM(WT, gKO, diag_net = 1, wij = "bij", symme = "FALSE", rm0 = "FALSE")
D <- test_fun(H0r, Hggm[[1]])
print(D[order(D$A), ])

