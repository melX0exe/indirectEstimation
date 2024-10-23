# Regression estimator

#' @param N Population size
#' @param x Vector of data
#' @param y Vector of auxiliary information
#' @param mY Mean of auxiliary variable
#' @export

Est.reg <- function(x,y,N,mY){
  
  n=length(x) # sample size
  s2x=var(x) # sample variance (x)
  s2y=var(y) # sample variance (y)
  sxy=cov(x,y) # sample covariance
  rho <- sxy/sqrt(s2x*s2y)
  
  # Estimators
  b <- sxy/s2y # Coefficient
  e.m <- mean(x) + b * (mY - mean(y)) # Mean
  e.t <- e.m * N # Total
  
  # Variance of the estimators
  v.m <- ((N-n)*s2x*(1-rho^2))/N*n # Mean
  v.t <- N^2*v.m # Total
  
  list(Est.b=b,
       Est.media=e.m,
       Est.total=e.t,
       Var.media=v.m,
       Var.total=v.t
       )
}