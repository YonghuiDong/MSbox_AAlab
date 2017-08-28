#' @title tanimoto similarity
#' @description calculate tanimoto similarities of different compounds
#' @param x similary
#' @export
#' @examples
#'  x <- data.frame(Samp1=c(0,0,0,1,1,1,0,0,1), Samp2=c(1,1,1,1,1,1,0,0,1))
#'  tanimoto(x)

tanimoto <- function(x) {
  res<-sapply(x, function(x1){
    sapply(x, function(x2) {i=length(which(x1 & x2)) / length(which(x1 | x2)); ifelse(is.na(i), 0, i)})
  })
  return(res)
}


