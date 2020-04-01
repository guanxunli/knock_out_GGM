## Daniel's method
pcCoefficients <- function(K, nComp) {
  # Taking out the gene to be regressed out
  y <- X[, K]
  Xi <- X
  Xi <- Xi[, -K]
  # Step 1: Perform PCA on the observed covariates data matrix to obtain $n$ number of the principal components.
  coeff <- RSpectra::svds(Xi, nComp)$v
  score <- Xi %*% coeff
  # Step 2: Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients.
  score <-
    Matrix::t(Matrix::t(score) / (apply(score, 2, function(X) {
      sqrt(sum(X ^ 2))
    }) ^ 2))
  # Step 3: Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
  Beta <- colSums(y * score)
  Beta <- coeff %*% (Beta)
  
  return(Beta)
}

## My method
pcCoefficients_l <- function(K, nComp) {
  # Taking out the gene to be regressed out
  y <- X[, K]
  Xi <- X
  Xi <- Xi[, -K]
  
  Xi_svd <- RSpectra::svds(Xi, nComp)
  
  return(t(t(Xi_svd$v)/Xi_svd$d) %*% crossprod(Xi_svd$u, y))
}

myid <- function(X, Y){
  return(sum((X - Y)^2))
}

X <- matrix(rnorm(1000000), nrow = 1000)
beta_d <- pcCoefficients(1, 3)
beta_l <- pcCoefficients_l(1, 3)
myid(beta_d, beta_l)
library(microbenchmark)
microbenchmark(
  pcCoefficients(2, 4),
  pcCoefficients_l(2, 4),
  times = 10
)