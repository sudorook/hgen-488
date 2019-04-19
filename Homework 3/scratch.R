#! /usr/bin/env Rscript

v <- c(10, 4, 2, 5, 6, 10)
w <- c(3, 2, 6, 10, 20, 4)

W <- 15

traceback <- function(P, w) {
  i = nrow(P)
  j = ncol(P)

  items <- c()
  while( (j > 0) && (i > 0) ) {
    if (P[i,j] == 1) {
      items <- append(items, i)
      j = j - w[i]
      i = i - 1
    } else {
      i = i - 1
    }
  }
  return(items)
}

knapsackSearch <- function(v, w, W) {
  n <- length(v)
  m <- matrix(0, ncol=(W+1), nrow=(n+1))
  p <- matrix(0, ncol=(W+1), nrow=(n+1))

  for (i in 2:(n+1)) {
    for (j in 2:(W+1)) {
      if (w[i-1] > j-1) {
        m[i,j] <- m[i-1,j]
        p[i,j] <- 0
      } else {
        m[i,j] <- max(m[i-1,j], m[i-1, j-w[i-1]] + v[i-1])
        if (m[i-1,j] < (m[i-1, j-w[i-1]] + v[i-1]))
          p[i,j] <- 1
        else
          p[i,j] <- 0
      }
    }
  }
  
  print((m[2:(n+1),2:(W+1)]))
  print((p[2:(n+1),2:(W+1)]))

  items <- traceback((p[2:(n+1),2:(W+1)]), w)
  res <- list(score=m[n+1,W+1], items=items)
  return(res)
}
