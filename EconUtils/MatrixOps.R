# MatrixOps.R
# Description: Some useful matrix operations, decompositions, and such.
# Author: Ari Boyarsky (aboyarsky@uchicago.edu)

# Calculate SVD
# Todo: Write my own eigenvalue decomp
my.svd <- function(dta){
  A <- as.matrix(dta)
  AtA <- t(A)%*%A
  AAt <- A%*%t(A)
  
  # find eigenvalues of AAt using eigen
  evd1 <- eigen(AAt)
  U <- evd1$vectors
  
  evd2 <- eigen(AtA)
  V <- evd2$vectors
  
  S <- diag(sqrt(evd2$values))
  
  return(list("U"= U, "V"= V, "S"= S))
}

# Demean using Mean Matrix
mat.demean <- function(dta){
  D <- as.matrix(dta)
  n <- NROW(D)
  O <- matrix(1, nrow=n, ncol=n)
  C <- diag(n) - (1 / n) * O
  return(as.matrix(C)%*%as.matrix(D))
}

# Angle between vectors
vector_angle <- function(a,b){
  theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
  return(theta)
}

# Find inverse of square matrix (must be positive definite -- admit to EVD)
inv.sq <- function(A){
  evd <- eigen(A)
  D <- diag(sqrt(evd$values))
  Q <- evd$vectors
  B <- solve(Q%*%D%*%t(Q))
  return(B)
}

# find X bar (mean) matrix
mean.mat <- function(X){
  X <- as.matrix(X)
  X_bar <- (1/NROW(X)) * t(X) %*% matrix(1, nrow = NROW(X), ncol = NCOL(X))
  
  return(X_bar)
  
}