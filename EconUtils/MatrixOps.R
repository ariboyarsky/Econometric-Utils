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
