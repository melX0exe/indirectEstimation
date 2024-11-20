# Ratio estimator

#' @param N Population size
#' @param x Vector of data
#' @param y Vector of auxiliary information
#' @param mY Mean of auxiliary variable
#' @param n1 Size of the sample (phase 1)
#' @param conf.level Confidence level
#' @param twophase Indicates whether double sampling was used
#' @export

Est.r<-function(x, y, N, mY, n1, conf.level=0.95, twophase=FALSE){
  
  n=length(x) # sample size
  s2x=var(x) # sample variance (x)
  s2y=var(y) # sample variance (y)
  sxy=cov(x,y) # sample covariance
  mx = mean(x)
  my = mean(y)
  
  # Estimators
  e.r=mx/my # Ratio
  e.m=mY*e.r # Mean
  e.t=e.m*N # Total
  
  # Variance of the estimators
  
  if (twophase == FALSE){ # We know the mean of the population for Y
    
    # Variance
    v.r=(N-n)*(s2x+e.r^2*s2y-2*e.r*sxy)/(N*n*mY^2) # Ratio
    v.m=v.r*mY^2 # Mean
    v.t=v.m*N^2 # Total
    
    # Errors
    er.r <- sqrt(v.r)
    er.m <- sqrt(v.m)
    er.t <- sqrt(v.t)
    
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
    
    return(list(Est.razon=e.r,
         Est.m=e.m,
         Est.t=e.t,
         ECM.razon=ecm.r,
         ECM.m=ecm.m,
         ECM.t=ecm.t,
         ratio.conf.int = c(lower.r,upper.r),
         mean.conf.int = c(lower.m,upper.m),
         total.conf.int = c(lower.t,upper.t)
    ))
    
  } else { # We used two-phased sampling
    
    # Variance
    v.m=(N-n1)*s2x/(N*n1) + (n1-n)/(n1*n)*(s2x+s2y*e.r^2-2*e.r*sxy) # Mean
    v.t=v.m*N^2 # Total
    
    # Errors
    er.m <- sqrt(v.m)
    er.t <- sqrt(v.t)
    
    # Confidence intervals
    # Mean
    lower.m <- e.m - qnorm(1-(1-conf.level)/2) * er.m
    upper.m <- e.m + qnorm(1-(1-conf.level)/2) * er.m
    # Total
    lower.t <- e.t - qnorm(1-(1-conf.level)/2) * er.t
    upper.t <- e.t + qnorm(1-(1-conf.level)/2) * er.t
    
    return(list(Est.razon=e.r,
         Est.m=e.m,
         Est.t=e.t,
         Var.m=v.m,
         Var.t=v.t,
         Error.m=er.m,
         Error.t=er.t,
         mean.conf.int = c(lower.m,upper.m),
         total.conf.int = c(lower.t,upper.t)
    ))
    
  }
  
}
