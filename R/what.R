#' @title search for m/z in from the idiom metabolomics database
#' @description tentative metabolite identification based on m/z value search
#' @author Yonghui Dong
#' @param mz  m/z value
#' @param ppm mass tolerance, default value = 10
#' @param mode ionization mode, either positive '+' or negative '-'
#' @export
#' @examples
#'  what(133.014, ppm = 10, mode = '-')

what <- function (mz, mode = c('+', '-'), ppm = 10) {

  ##(1) input check
  if(is.numeric(mz) == FALSE) {stop("warning: mass to charge ratio mz shoule be numeric!")}
  if(mode != "+" & mode !="-") {stop("warning: ion mode invalid. Choose '+' or '-'.\n")}

  ##(2) load database
  DB <- as.data.frame(sysdata$plant_db)

  ##(3) search in database
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))
  Result <- vector("list", length(mz))

  if(mode == '-') {
    for (i in 1:length(mz)) {
      DB.list <- expand.grid.df(mz[i], DB[, -(3:5)])
      colnames(DB.list)[1] <- "search"
      cal_ppm <- with(DB.list, (DB.list$`[M-H]` - search) * 10^6 / DB.list$`[M-H]`)
      cal_ppm <-  round(cal_ppm, digits = 2)
      DB.list <- cbind(DB.list, ppm = cal_ppm)
      Result[[i]] = DB.list[(abs(cal_ppm) <= ppm), ]
    }
  } else {
    for (i in 1:length(mz)) {
      ## get [M+H]
      DB.list_H <- expand.grid.df(mz[i], "[M+H]+", DB[, -(4:6)])
      colnames(DB.list_H)[c(1, 2, 5)] <- c("search", "Adduct", "mzs")
      ## get [M+Na]
      DB.list_Na <- expand.grid.df(mz[i], "[M+Na]+", DB[, -c(3, 5, 6)])
      colnames(DB.list_Na)[c(1, 2, 5)] <- c("search", "Adduct", "mzs")
      ## get [M+K]
      DB.list_K <- expand.grid.df(mz[i], "[M+K]+", DB[, -c(3, 4, 6)])
      colnames(DB.list_K)[c(1, 2, 5)] <- c("search", "Adduct", "mzs")
      ## combine them
      DB.list <- rbind(DB.list_H, DB.list_Na, DB.list_K)
      cal_ppm <- with(DB.list, (mzs - search) * 10^6 / mzs)
      cal_ppm <-  round(cal_ppm, digits = 2)
      DB.list <- cbind(DB.list, ppm = cal_ppm)
      Result[[i]] = DB.list[(abs(cal_ppm) <= ppm), ]
    }
  }
  search_result <- do.call(rbind.data.frame, Result)

  ##(4) check if Result is empty, and return result
  if(nrow(search_result) == 0) {
    message('Not Found, Unknown')
  } else {
    return(search_result)
  }
}
