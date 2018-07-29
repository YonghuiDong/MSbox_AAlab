#' @title search for m/z in home database
#' @description tentative metabolite identification based on m/z value search in a home database
#' @author Yonghui Dong
#' @param mz  m/z value
#' @param ppm mass tolerance, default value = 10
#' @param mode ionization mode, either positive '+' or negative '-'
#' @param relatives, shuould the corresponding adducts, fragments and parent ions shown
#' @export
#' @examples
#'  what(133.014, ppm = 10, mode = '-')
#'  what(133.014, ppm = 10, mode = '-', relative = TRUE)

what <- function (mz, mode = c('+', '-'), ppm = 10, relatives = FALSE) {

  ##(1) input check
  if(is.numeric(mz) == FALSE) {stop("Warning: mass to charge ratio mz shoule be numeric!")}
  if(mode != "+" & mode !="-") {stop("WARNING: ion mode invalid. Choose '+' or '-'.\n")}

  ##(2) load database
  DB_pos <- as.data.frame(sysdata$analyticon_pos)
  DB_neg <- as.data.frame(sysdata$analyticon_neg)

  ##(3) search in database
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL),
                                         list(...))
  if(mode == '-') {
    DB.list <- expand.grid.df(mz, DB_neg)
    colnames(DB.list)[1] <- 'my_mz'
    myppm <- with(DB.list, abs(my_mz - mz) * 10^6 / mz)
    Result = DB.list[(myppm <= ppm), ]
  } else {
    DB.list <- expand.grid.df(mz, DB_pos)
    colnames(DB.list)[1] <- 'my_mz'
    myppm <- with(DB.list, abs(my_mz - mz) * 10^6 / mz)
    Result = DB.list[(myppm <= ppm), ]
  }

  ##(4) check if Result is empty
  if(nrow(Result) == 0) stop('Not Found, Unknown')

  ##(5) return result
  if(nrow(Result) != 0) {
   ## include rows which have the same ID as the selected m/z (same parent ion).
    search_result <- DB.list[DB.list$ID %in% Result$ID,]
  } else
    message('Not Found, unknown')
  ifelse(relatives == TRUE, return(search_result), return(Result))
}
