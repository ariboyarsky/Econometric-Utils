# UnsupervisedLearning.R
#
# Description: Common analysis techniques often categorized as 
# unsupervised learning.
#
# Author: Ari Boyarsky (aboyarsky@uchicago.edu)

# Set working dir to current dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load Matrix Ops Functions
source("MatrixOps.R")

#---------
# Factor Analysis
#---------
factor_analysis <- function(M, type="Cov", dta){
  # Function preforms factor analysis finding the 
  # loading matrix. 
  # Arguments: 
  #   M := m<p number otf variables from data
  #   type := "Cov" or "Cor" matrix method
  #   dta := dataframe
  
  # E = S - LL^T - psi
  if(type == "Cov"){
    S <- 1/(NROW(dta)-1) * t(demean(dta))%*%demean(dta)
  }else{
    S <- cor(dta)
  }
  # EVD of S
  evd <- eigen(S)
  Q <- evd$vectors[,1:M] # is this the correct dim?
  Lambda <- evd$values[1:M]
  
  # Calculate L = Q*Lambda^{1/2}
  L1 <- sqrt(Lambda)*Q
  # Calc psi diag(S-LL^T)
  psi <- diag(S - L1%*%t(L1))
  
  variability <- evd$values / sum(S)
  
  return(list("L" = L1, 
              "psi" = psi,
              "var" = variability))
}

#---------
# Procrustes Analysis
#---------

procrustes <- function(A,B){
  # Function preforms procrustes analysis.
  # Find Q matrix that rotates B s.t. Frobeneous norm is minimized.
  # Arguments:
  #   A := mxn matrix
  #   B := mxn matrix
  
  # procrustes analysis := min_Q ||A-BQ||_F == max tr(Q^Tb^TA)
  # Q = UV^T where C = B^TA = UDV^T
  C <- t(B)%*%A
  
  # SVD of C
  s <- svd(C)
  
  # right/left sing vectors
  u <- s$u
  v <- s$v
  
  # return Q = UV^T
  return(u%*%t(v))
}

#---------
# Canonical Correlation Analysis
#---------

smpl.cancor <- function(X,Y, type = "cov"){
  if(type == "cov"){
    # sample covariance version
    Sx <- 1/(NROW(X)-1)*t(demean(X))%*%demean(X)
    Sy <- 1/(NROW(Y)-1)*t(demean(Y))%*%demean(Y)
    Sxy <- 1/(NROW(Y)-1)*t(demean(X))%*%demean(Y)
  }else{
    # sample correlation version
    Sx <- cor(X)
    Sy <- cor(Y)
    Sxy <- cor(X,Y)
  }
  
  
  Sx.inv <- inv.sq(Sx)
  Sy.inv <- inv.sq(Sy)
  
  G <- Sx.inv %*% Sxy %*% Sy.inv
  decomp <- svd(G)
  
  r <- decomp$d
  
  a <- Sx.inv %*% decomp$u
  b <- Sy.inv %*% decomp$v
  
  return(list(
    "cor" = r,
    "a" = a,
    "b" = b
  ))
}