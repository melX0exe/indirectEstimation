# Compare Ratio and Regression estimators

#' @param x Vector of data
#' @param y Vector of auxiliary information
#' @param b Slope of the regression estimator
#' @export

compare.reg <- function(x,y,b=NA){
  
  if(is.na(b)){
    b <- cov(x,y)/var(y) # Slope
  }
  
  r <- mean(x)/mean(y) # Ratio
  
  return(ifelse(all.equal(r,b)==TRUE,"Equivalent", "Regression"))
  
}