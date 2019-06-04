set.seed(10)

# Biclustering demonstration
library(biclust)

# n <- 20
# m <- 20

# mus <- c(-1, 0, 1)
# sds <- c(.1, .1, .1)
# 
# X <- matrix(0, nrow=m, ncol=n)
# for (i in 1:m) {
#   for (j in 1:n) {
#     id <- ceiling(runif(1)*4)
#     if ( id <= length(mus) ) {
#       X[i,j] <- rnorm(1, mus[id], sds[id])
#     }
#   }
# }

# X <- matrix(rbinom(n*m, 50, c(.6, .4, .2)), nrow=m, ncol=n)
# # res <- biclust(X, method=BCXmotifs(), alpha=0.05, number=50)
# res <- biclust(X, method=BCCC(), alpha=0.05, number=50)
# 
# heatmap(x=X, Rowv=NA, Colv=NA, reorderfun=NA)
# heatmapBC(x=X, bicResult=res, xlab=NA, ylab=NA)


data(BicatYeast)
X <- discretize(BicatYeast)
# X <- BicatYeast

X.df <- data.frame(X)
res <- biclust(X, method=BCXmotifs(), alpha=0.05, number=50)
# res <- biclust(X, method=BCCC(), delta=1.5, alpha=0.05, number=50)
res <- biclust(X, method=BCCC(), alpha=0.05, number=10)
# res2 <- biclust(X, method=BCCC())

png('heatmap.png')
heatmap(x=X, Rowv=NA, Colv=NA, reorderfun=NA)
dev.off()

png('heatmapBC.png')
heatmapBC(x=X, bicResult=res, xlab=NA, ylab=NA)
dev.off()


# Sparse SVD

library('ssvd')

# X <- matrix(rnorm(2^15), 2^7, 2^8)
X <- matrix(rnorm(100), 10, 10)
X[sample(1:length(X), 10)] <- rnorm(10, 5, 2)
res <- ssvd(X, method = "method")
ans.initial <- ssvd.initial(X, method = "method")
ans.iter <- ssvd.iter.thresh(X, u.old=ans.initial$u, v.old=ans.initial$v, method = "method")

#

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
