#' @title search for m/z in from the idiom metabolomics database
#' @description tentative metabolite identification based on m/z value search
#' @author Yonghui Dong
#' @param mz  m/z values
#' @param ppm mass tolerance, default value = 10
#' @param mode ionization mode, either positive '+' or negative '-'
#' @export
#' @examples
#'  what(133.014, ppm = 10, mode = '-')

what <- function (mz, mode = NULL, ppm = 5) {

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
      width <- options()$width * 0.3
      cat(paste0(rep(c(intToUtf8(0x2698), "="), i / length(mz) * width), collapse = ''))
      cat(paste0(round(i / length(mz) * 100), '% completed'))
      DB.list <- expand.grid.df(i, mz[i], DB[, -(3:5)])
      colnames(DB.list)[c(1, 2)] <- c("QueryID", "Search")
      cal_ppm <- with(DB.list, (DB.list$`[M-H]` - Search) * 10^6 / DB.list$`[M-H]`)
      cal_ppm <-  round(cal_ppm, digits = 2)
      DB.list <- cbind(DB.list, ppm = cal_ppm)
      Result[[i]] = DB.list[(abs(cal_ppm) <= ppm), ]
      row.names(Result[[i]]) <- NULL
      if (i == length(mz)) cat(': Searching Done.')
      else cat('\014')
    }
  } else {
    for (i in 1:length(mz)) {
      width <- options()$width * 0.3
      cat(paste0(rep(c(intToUtf8(0x2698), "="), i / length(mz) * width), collapse = ''))
      cat(paste0(round(i / length(mz) * 100), '% completed'))
      ## get [M+H]
      DB.list_H <- expand.grid.df(i, mz[i], "[M+H]+", DB[, -(4:6)])
      colnames(DB.list_H)[c(1, 2, 3, 6)] <- c("QueryID", "Search", "Adduct", "mzs")
      ## get [M+Na]
      DB.list_Na <- expand.grid.df(i, mz[i], "[M+Na]+", DB[, -c(3, 5, 6)])
      colnames(DB.list_Na)[c(1, 2, 3, 6)] <- c("QueryID", "Search", "Adduct", "mzs")
      ## get [M+K]
      DB.list_K <- expand.grid.df(i, mz[i], "[M+K]+", DB[, -c(3, 4, 6)])
      colnames(DB.list_K)[c(1, 2, 3, 6)] <- c("QueryID", "Search", "Adduct", "mzs")
      ## combine them
      DB.list <- rbind(DB.list_H, DB.list_Na, DB.list_K)
      cal_ppm <- with(DB.list, (mzs - Search) * 10^6 / mzs)
      cal_ppm <-  round(cal_ppm, digits = 2)
      DB.list <- cbind(DB.list, ppm = cal_ppm)
      Result[[i]] = DB.list[(abs(cal_ppm) <= ppm), ]
      row.names(Result[[i]]) <- NULL
      if (i == length(mz)) cat(': Searching Done.')
      else cat('\014')
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
