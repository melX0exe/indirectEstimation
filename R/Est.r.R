# Ratio estimator

#' @param N Population size
#' @param x Vector of data
#' @param y Vector of auxiliary information
#' @param mY Mean of auxiliary variable
#' @param conf.level Confidence level
#' @export

Est.r<-function(x, y, N, mY, conf.level=0.95){
  
  n=length(x) # sample size
  s2x=var(x) # sample variance (x)
  s2y=var(y) # sample variance (y)
  sxy=cov(x,y) # sample covariance
  
  # Estimators
  e.r=mean(x)/mean(y) # Ratio
  e.m=mY*e.r # Mean
  e.t=e.m*N # Total
  
  # Variance of the estimators
  v.r=(N-n)*(s2x+e.r^2*s2y-2*e.r*sxy)/(N*n*mY^2) # Ratio
  v.m=v.r*mY^2 # Mean
  v.t=v.m*N^2 # Total
  
  # Errors
  er.r <- sqrt(v.r)
  er.m <- sqrt(v.m)
  er.t <- sqrt(v.m)
  
  # Bias
  b.r=(N-n)*(e.r*s2y-sxy)/(N*n*mY^2) # Ratio
  b.m=b.r*mY # Mean
  b.t=b.m*N # Total
  
  # MSE 
  ecm.r=v.r+b.r^2 # Ratio
  ecm.m=v.m+b.m^2 # Mean
  ecm.t=v.t+b.t^2 # Total
  
  # Confidence intervals
  # Ratio
  lower.r <- e.r - b.r - qnorm(1-(1-conf.level)/2) * er.r
  upper.r <- e.r - b.r + qnorm(1-(1-conf.level)/2) * er.r
  # Mean
  lower.m <- e.m - b.m - qnorm(1-(1-conf.level)/2) * er.m
  upper.m <- e.m - b.m + qnorm(1-(1-conf.level)/2) * er.m
  # Total
  lower.t <- e.t - b.t - qnorm(1-(1-conf.level)/2) * er.t
  upper.t <- e.t - b.t + qnorm(1-(1-conf.level)/2) * er.t
  
  list(Est.razon=e.r,
       Est.m=e.m,
       Est.t=e.t,
       ECM.razon=ecm.r,
       ECM.m=ecm.m,
       ECM.t=ecm.t,
       ratio.conf.int = c(lower.r,upper.r),
       mean.conf.int = c(lower.m,upper.m),
       total.conf.int = c(lower.t,upper.t)
       )
}
