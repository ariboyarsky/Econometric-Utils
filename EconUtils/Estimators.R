# Estimators.R
# Description: Functions for estimating parameters in a variety of 
#                 economic structural and reduced form models.
# Author: Ari Boyarsky (aboyarsky@uchicago.edu)

##############################################################################
# logit regression with graident optim method 
# using linear regression to get init values as a nice trick
##############################################################################

log.logit <- function(X,Y){
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
  
  return(logitmodel)
}
##############################################################################
##############################################################################

# Mutliple Linear Regression
# Homoskedastic Std Errors (todo: add option for asymptotically consistent SE)
linreg <- function(Y, X, con = 1){
  # get stuff as matricies
  Y <- as.matrix(Y)
  
  if(con == 1){
    X <- as.matrix(cbind(1,X))
  }else{
    X <- as.matrix(X)
  }

  
  # helpers
  XtX <- t(X)%*%X
  XtY <- t(X)%*%Y
  
  beta <- solve(XtX)%*%XtY
  
  # get SE
  pred <- X%*%beta
  n <- NROW(X)
  sigma2 <- sum((Y-pred)^2) / (nrow(X)-ncol(X))
  SE <- sqrt(diag(solve(XtX)*sigma2))
  
  return(list("beta"=beta,"SE"=SE))
}
##############################################################################
##############################################################################
# IV Regression (2SLS) 
# Consider adding asymptotically consistent SE
# ^Doing so yields a memory issue in R with the Angrist data???
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
  
  # se
  p2 <- pred%*%beta.iv
  n <- NROW(X)
  sigma2 <- sum((Y-p2)^2) / (nrow(pred)-ncol(pred))
  SE <- sqrt(diag(solve(t(pred)%*%pred)*sigma2))
  
  return(list("beta.iv"=beta.iv, "beta.first"=beta.first, "SE"=SE))
}
##############################################################################
##############################################################################

