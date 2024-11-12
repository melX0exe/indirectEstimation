# Simple random sampling estimator

#' @param N Population size
#' @param x Vector of data
#' @param conf.level Confidence level
#' @export


Est.mas<-function(x, N, conf.level=0.95){
  
  n=length(x) # sample size
  s2=var(x) # sample variance
  
  # Estimators
  e.m=mean(x) # Mean estimator
  e.t=e.m*N # Total estimator
  
  # Variance of the estimators
  v.m=(N-n)*s2/(N*n) # Mean
  v.t=v.m*N^2 # Total
  
  # Errors
  er.m=sqrt(v.m) # Mean
  er.t=sqrt(v.t) # Total
  
  # Confidence intervals
  # Mean
  lower.m <- e.m - qnorm(1-(1-conf.level)/2) * er.m
  upper.m <- e.m + qnorm(1-(1-conf.level)/2) * er.m
  # Total
  lower.t <- e.t - qnorm(1-(1-conf.level)/2) * er.t
  upper.t <- e.t + qnorm(1-(1-conf.level)/2) * er.t
  
  # Output
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
