# Regression estimator

#' @param N Population size
#' @param x Vector of data
#' @param y Vector of auxiliary information
#' @param mY Mean of auxiliary variable
#' @param n1 Size of the sample (phase 1)
#' @param b Slope
#' @param conf.level Confidence level
#' @param twophase Indicates whether double sampling was used
#' @export

Est.reg <- function(x, y, N, mY, n1, b=NA, conf.level=0.95, twophase=FALSE){
  
  n=length(x) # sample size
  s2x=var(x) # sample variance (x)
  s2y=var(y) # sample variance (y)
  sxy=cov(x,y) # sample covariance
  rho <- sxy/sqrt(s2x*s2y)
  
  # Estimators
  if(is.na(b)){
    b <- sxy/s2y # Coefficient
  }
  e.m <- mean(x) + b * (mY - mean(y)) # Mean
  e.t <- e.m * N # Total
  
  if(twophase == FALSE){ # We know the mean of the population for Y
    
    # Variance of the estimators
    v.m <- ((N-n)*s2x*(1-rho^2))/(N*n) # Mean
    v.t <- N^2*v.m # Total
    
    # Errors
    er.m=sqrt(v.m) # Mean
    er.t=sqrt(v.t) # Total
    
    # Confidence interval
    # We assume bias = 0:
    # Mean
    lower.m <- e.m - qnorm(1-(1-conf.level)/2) * er.m
    upper.m <- e.m + qnorm(1-(1-conf.level)/2) * er.m
    # Total
    lower.t <- e.t - qnorm(1-(1-conf.level)/2) * er.t
    upper.t <- e.t + qnorm(1-(1-conf.level)/2) * er.t
    
    return(list(Est.b=b,
         Est.media=e.m,
         Est.total=e.t,
         Var.media=v.m,
         Var.total=v.t,
         Err.m=er.m,
         Err.t=er.t,
         mean.conf.int = c(lower.m,upper.m),
         total.conf.int = c(lower.t,upper.t)
    ))
    
  } else {
    
    # Variance of the estimators
    v.m <- (N-n1)*s2x/(N*n1) + (n1-n)/(n1*n)*(1-rho^2)*s2x # Mean
    v.t <- N^2*v.m # Total
    
    # Errors
    er.m=sqrt(v.m) # Mean
    er.t=sqrt(v.t) # Total
    
    # Confidence interval
    # We assume bias = 0:
    # Mean
    lower.m <- e.m - qnorm(1-(1-conf.level)/2) * er.m
    upper.m <- e.m + qnorm(1-(1-conf.level)/2) * er.m
    # Total
    lower.t <- e.t - qnorm(1-(1-conf.level)/2) * er.t
    upper.t <- e.t + qnorm(1-(1-conf.level)/2) * er.t
    
    return(list(Est.b=b,
         Est.media=e.m,
         Est.total=e.t,
         Var.media=v.m,
         Var.total=v.t,
         Err.m=er.m,
         Err.t=er.t,
         mean.conf.int = c(lower.m,upper.m),
         total.conf.int = c(lower.t,upper.t)
    ))
    
  }
  
}
