#' @title Contaminants in MS
#' @description check the possible contaminants
#' @param mz suspected m/z value
#' @param ppm mass tolerance, default value = 10
#' @param mode ionization mode, either positive '+' or negative '-'
#' @export
#' @examples
#'  contam(33.0335, ppm = 10, mode = '+')
#'  contam(44.998, ppm = 10, mode = '-')

contam <- function (mz, ppm = 10, mode = c('+', '-')) {
  if(is.numeric(mz) == FALSE) {stop("Warning: mass to charge ratio mz shoule be numeric!")}
  if(mode != "+" & mode !="-") {stop("WARNING: ion mode invalid. Choose '+' or '-'.\n")}
  contam_pos <- as.data.frame(sysdata$contam_pos)
  contam_neg <- as.data.frame(sysdata$contam_neg)
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL),
                                         list(...))
  if(mode == '-') {
    contam.list <- expand.grid.df(mz, contam_neg)
    colnames(contam.list)[1] <- 'mymz'
    myppm <- with(contam.list, abs(mymz - mz) * 10^6 / mz)
    Result = contam.list[(myppm <= ppm), -1]
    if(nrow(Result) != 0) {
      return(Result)
    } else
      message('Not Found')
  } else if (mode == '+') {
    contam.list <- expand.grid.df(mz, contam_pos)
    colnames(contam.list)[1] <- 'mymz'
    myppm <- with(contam.list, abs(mymz - mz) * 10^6 / mz)
    Result = contam.list[(myppm <= ppm), -1]
    if(nrow(Result) != 0) {
      return(Result)
    } else
      message('Not Found, Unknown')
   }
 }


