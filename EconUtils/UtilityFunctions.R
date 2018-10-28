# UtilityFunctions.R
# Description: Small utility functions useful in statistical analysis.
# Author: Ari Boyarsky (aboyarsky@uchicago.edu)

# demean a dta
demean <- function(dta){
  demeaned <- dta
  for(i in 1:ncol(dta)){
    m <- mean(dta[,i])
    for(j in 1:nrow(dta)){
      demeaned[j,i] <- dta[j,i] - m
    }
  }
  return(demeaned)
}

# scale by sd
sd_scale <- function(dta){
  scaled <- dta
  for(i in 1:ncol(dta)){
    s <- sd(dta[,i])
    for(j in 1:nrow(dta)){
      scaled[j,i] <- dta[j,i]/s
    }
  }
  return(scaled)
}