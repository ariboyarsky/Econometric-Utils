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