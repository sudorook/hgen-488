#! /usr/bin/env Rscript

v <- c(0, 10, 4, 2, 5, 6, 10)
w <- c(0, 3, 2, 6, 10, 20, 4)

W <- sum(w)
n <- length(v)

m <- matrix(0, ncol=(W+1), nrow=n)
p <- matrix(0, ncol=(W+1), nrow=n)

for (i in 2:(n)) {
  for (j in 1:(W+1)) {
    if (w[i] > j) {
      m[i,j] <- m[i-1,j]
      p[i,j] <- -1
    } else {
      m[i,j] <- max(m[i-1,j], m[i-1, j-w[i]] + v[i])
      if (length(m[i-1, j-w[i]] + v[i]) == 0) {
        p[i,j] <- -1
      } else if (m[i-1,j] < (m[i-1, j-w[i]] + v[i])) {
        p[i,j] <- 1
      } else {
        p[i,j] <- -1
      }
    }
  }
}

t(m[2:(n),2:(W+1)])
t(p[2:(n),2:(W+1)])
