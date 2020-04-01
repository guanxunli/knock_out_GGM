set.seed(1234)
X <- matrix(rnorm(100), nrow = 10)
K <- 1
nComp <- 3

y <- X[, K]
Xi <- X
Xi <- Xi[, -K]
p <- nrow(Xi)
q <- ncol(Xi)

# Xi_svd <- RSpectra::svds(Xi, nComp)
Xi_svd <- svd(Xi)
u <- Xi_svd$u
v <- Xi_svd$v
d <- Xi_svd$d

a <- X[, 1] - X[, 2]
b <- c(1, rep(0, q - 1))

m <- crossprod(u, a)
p <- a - u %*% m
Ra <- sqrt(sum(p^2))
P <- p/Ra

n <- crossprod(v, b)
q <- b - v %*% n
Rb <- sqrt(sum(q^2))
Q <- q/Rb

K <- diag(c(d, 0)) + tcrossprod(c(m, Ra), c(n, Rb))

sol1 <- cbind(v[, 1:3], q) %*% solve(K[1:4, 1:4], t(cbind(u[, 1:3], p)))

K <- 2
y <- X[, K]
Xi <- X
Xi <- Xi[, -K]
p <- nrow(Xi)
q <- ncol(Xi)

Xi_svd <- RSpectra::svds(Xi, nComp)
# Xi_svd <- svd(Xi)
u <- Xi_svd$u
v <- Xi_svd$v
d <- Xi_svd$d

sol2 <- t(t(Xi_svd$v)/Xi_svd$d) %*% t(u)

myid <- function(X, Y){
  return(sum((X - Y)^2))
}
myid(sol1, sol2)
