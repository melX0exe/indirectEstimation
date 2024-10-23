# Difference estimator

#' @param N Population size
#' @param x Vector of data
#' @param y Vector of auxiliary information
#' @param mY Mean of auxiliary variable
#' @param conf.level Confidence level
#' @export

Est.dif<-function(x,y,N,mY,conf.level=0.95){
  
  n=length(x) # sample size
  s2x=var(x) # sample variance (x)
  s2y=var(y) # sample variance (y)
  sxy=cov(x,y) # sample covariance
  
  # Estimators
  e.m=mY+mean(x)-mean(y) # Mean
  e.t=e.m*N # Total
  
  # Variance of the estimators
  v.m=(N-n)*(s2x+s2y-2*sxy)/(N*n) # Mean
  v.t=v.m*N^2 # Total
  
  # Errors
  er.m=sqrt(v.m) # Mean
  er.t=sqrt(v.t) # Sample
  
  # Confidence intervals
  # Mean
  lower.m <- e.m - qnorm(1-(1-conf.level)/2) * er.m
  upper.m <- e.m + qnorm(1-(1-conf.level)/2) * er.m
  # Total
  lower.t <- e.t - qnorm(1-(1-conf.level)/2) * er.t
  upper.t <- e.t + qnorm(1-(1-conf.level)/2) * er.t
  
  list(Est.m=e.m,
       Est.t=e.t,
       Var.m=v.m,
       Var.t=v.t,
       Err.m=er.m,
       Err.t=er.t,
       mean.conf.int = c(lower.m,upper.m),
       total.conf.int = c(lower.t,upper.t)
  )
}