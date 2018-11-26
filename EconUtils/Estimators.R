# Estimators.R
# Description: Functions for estimating parameters in a variety of 
#                 economic structural and reduced form models.
# Author: Ari Boyarsky (aboyarsky@uchicago.edu)

##############################################################################
# logit regression with graident optim method 
# using linear regression to get init values as a nice trick
##############################################################################

log.logit <- function(X,Y){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  # drop all 0 cols in X
  X <-X[, colSums(X != 0) > 0]
  
  # get start values by ols to speed up process
  startvalues <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  # logit loglikli
  K <- as.numeric(ncol(X))
  logit.logl <- function(beta) {
    
    # logit
    prob <- exp(X%*%beta)/(1+exp(X%*%beta)) 
    logprob <- log(prob)
    
    y0 <- 1 - Y
    logprob0 <- log(1 - prob)
    
    # get loglikli
    loglikli <- -sum(t(Y)%*%logprob + t(y0)%*%logprob0)
    return(loglikli)
  }
  
  # logit gradient function
  logit.grad <- function(beta) {
    grad <- beta*0
    
    # gradient
    prob <- exp(X%*%beta) /(1+exp(X%*%beta)) 
    
    # for each var get grad
    for (k in 1:K) { 
      grad[k] <- sum(X[,k] * (Y - prob))
    }
    return(-grad)
  }
  
  # Using BFGS method --- reqs report = 1
  logitmodel <- optim(startvalues, logit.logl, gr=logit.grad,
                      control=list(trace=FALSE, REPORT=1), 
                      method="BFGS", hessian=TRUE)
  
  # return pred values too
  logit.coeffs <- logitmodel$par
  
  # get predicted values
  pred <- X%*%logit.coeffs
  P <- exp(pred)/(1+exp(pred))
  
  return(list("Model"=logitmodel, "Pred"=P))
}
##############################################################################
##############################################################################

# Mutliple Linear Regression
# Options for weighted lest squares, robust se, and constant
linreg <- function(Y, X, con = 1, weights = NULL, hetero = TRUE){
  # get stuff as matricies
  Y <- as.matrix(Y)
  
  if(con == 1){
    X <- as.matrix(cbind(1,X))
  }else{
    X <- as.matrix(X)
  }
  
  # If weights submitted then calculate weighted OLS
  if(!is.null(weights)){
    W <- unlist(c(weights))
  } else{
    W <- replicate(nrow(X), 1)
  }
  XW <- X * W
  
  # remove any 0 cols
  X <-X[, colSums(X != 0) > 0]
  XW <-XW[, colSums(XW != 0) > 0]
  
  beta <- solve(t(XW) %*% X)%*%t(XW)%*%Y
  
  # residuals
  pred <- X%*%beta
  n <- NROW(X)
  k <- ncol(X)
  U <- Y - pred
  
  # get SE
  if(hetero == FALSE){
    sigma2 <- sum((Y-pred)^2) / (n-k)
    SE <- sqrt(diag(solve(t(XW)%*%X)*sigma2))
  }else{
    Omega <- Diagonal(length(U) ,as.numeric(U ^ 2))
    V <-  solve(t(XW) %*% X) %*% t(XW) %*% Omega %*% X %*% solve(t(XW) %*% X)
    SE <- sqrt(diag(V))
  }
  
  return(list("beta"=beta,"SE"=SE,"res"=U))
}
##############################################################################
##############################################################################
# IV Regression (2SLS) 
ols.iv <- function(Y,Z,X,controls, con = 1){
  if (con==1){
    controls <- as.matrix(cbind(1, controls))
  }else{
    controls <- as.matrix(controls) 
  }
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  
  # first stage
  W <- as.matrix(cbind(controls,Z))
  
  # drop any 0 cols
  W <-W[, colSums(W != 0) > 0]
  
  beta.first <- solve(t(W)%*%W, tol = 1e-20)%*%t(W)%*%X
  pred <- W%*%beta.first
  
  colnames(pred) <- c("Z")
  
  # second stage
  if(con==1){
    pred <- as.matrix(cbind(1,pred,C))
  }else{
    pred <- as.matrix(cbind(pred,C))
  }
  
  # drop any 0 cols
  pred <-pred[, colSums(pred != 0) > 0]
  
  
  beta.iv <- solve(t(pred)%*%pred)%*%t(pred)%*%Y
  
  # robust se
  p2 <- pred%*%beta.iv
  U <- Y - p2
  Omega <- Diagonal(length(U) ,as.numeric(U ^ 2))
  V = solve(t(pred)%*%pred)%*%t(pred)%*%Omega%*%pred%*%solve(t(pred)%*%pred)
  SE <- sqrt(diag(V))
  
  return(list("beta.iv"=beta.iv, "beta.first"=beta.first, "SE"=SE))
}
##############################################################################
##############################################################################

