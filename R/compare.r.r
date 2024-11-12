# Compare Ratio and s.r.s estimators

#' @param x Vector of data
#' @param y Vector of auxiliary information
#' @export

compare.r <- function(x,y){
  return(ifelse(cor(x,y) > (sqrt(var(y))*mean(x))/(sqrt(var(x))*mean(y)*2),"Ratio","s.r.s"))
}

