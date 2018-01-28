#' Test association of variables potentially given others 
#' 
#' Test association of variables potentially given others and return a p-value
#' 
#' @param x A numeric vector
#' @param y A numeric vector
#' @param z A numeric vector or matrix
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less". 
#'   You can specify just the initial letter. "greater" corresponds to positive association, 
#'   "less" to negative association.
#' @return A p-value
#' @details z is optional. If given, partial correlation test \code{\link[ppcor]{pcor.test}} is used.   

assoc = function(x, y, z=NULL, alternative = c("two.sided", "less", "greater")){
  alternative = match.arg(arg=alternative, choices=c("two.sided", "less", "greater"))
  if( is.null(z) ){
    pv = cor.test(x,y,alternative = alternative)$p.value
  } else {
    pcor.v = ppcor::pcor.test(x,y,z=z)
    if (alternative=="two.sided"){
      pv = pcor.v$p.value
    } else {
      if ((alternative=="greater" && pcor.v$estimate > 0)||(alternative=="less" && pcor.v$estimate < 0)){
        pv = pcor.v$p.value/2
      } else {
        pv = 1 - pcor.v$p.value/2
      }
    }
  }#end else !is.null(z)
  return(pv)
}