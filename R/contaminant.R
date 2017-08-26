#' @title Contaminants in MS
#' @description Check the possible contaminants
#' @param mz suspected m/z value
#' @param ppm mass tolerance, default value = 10
#' @param mode ionization mode, either positive '+' or negative '-'
#' @examples
#' contam(33.033, ppm = 10, mode = '+')
#' contam(44.998, ppm = 10, mode = '-')

contam <- function (mz, ppm = 10, mode = c('+', '-')) {
  data('contam_pos')
  data('contam_neg')
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL),
                                         list(...))
  if(mode == '-') {
    contam.list <- expand.grid.df(mz, contam_neg)
    colnames(contam.list)[1] <- 'mymz'
    myppm <- with(contam.list, abs(mymz - mz) * 10^6 / mz)
    Result = contam.list[(myppm <= ppm), -1]
    return(Result)
  } else if (mode == '+') {
    contam.list <- expand.grid.df(mz, contam_pos)
    colnames(contam.list)[1] <- 'mymz'
    myppm <- with(contam.list, abs(mymz - mz) * 10^6 / mz)
    Result = contam.list[(myppm <= ppm), -1]
    return(Result)
  } else
    message ('Ion mode is not correct')
}


