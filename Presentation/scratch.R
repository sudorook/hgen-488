set.seed(10)

n <- 10
m <- 10
ct <- .25

X <- matrix(0, ncol=n, nrow=m)
X[sample(seq(1,m*n), ct*n*m)] <- rnorm(ct*n*m)

res <- svd(X)

U <- res$u
V <- res$v
S <- diag(res$d)

round(res$d[1] * U[,1] %*% t(V[,1]), 5)
round(res$d[2] * U[,2] %*% t(V[,2]), 5)
round(res$d[3] * U[,3] %*% t(V[,3]), 5)
round(res$d[4] * U[,4] %*% t(V[,4]), 5)
round(res$d[5] * U[,5] %*% t(V[,5]), 5)

X_hat <- matrix(0, ncol=n, nrow=m)

for(i in 1:length(res$d)) {
  X_hat <- X_hat + res$d[i] * U[,i] %*% t(V[,i])
}
