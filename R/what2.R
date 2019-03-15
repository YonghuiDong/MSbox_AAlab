#' @title metabolite/lipid identification
#' @description metabolite/lipid identification using home database accoring to m/z and RT. MSMS fragments are also given.
#' @author Yonghui Dong
#' @param xset xcms object
#' @param RT.DB correction RT result from RT_correct()
#' @param use.DB which database is used for peak identification, choose between "Lipidomics" and "Metabolomics".
#' @param RT.cor should use corrected RT
#' @param RT.cor.method which corrected RT should be used? L: linear regression (default) or P: polynomial regression?
#' @param frag should search fragment
#' @param ppm mass tolerance, default value = 10
#' @param rt retention time search window, default value = 25
#' @importFrom xcms peakTable
#' @export
#'@examples
#'\dontrun{
#'my_compound <- what2(xset)
#'}

what2 <- function (xset, RT.DB = NULL, use.DB = "METABOLOMICS", RT.cor = TRUE, RT.cor.method = "L", frag = T, ppm = 10, rt = 25) {

  #(1) input check
  if(class(xset) != "xcmsSet") {stop("the input object is not an xcmsSet object")}
  if(!is.numeric(ppm)){stop("invalid calss of ppm threshold: not numeric")}
  if(!is.numeric(rt)){stop("invalid calss of RT threshold: not numeric")}
  if(ppm < 0){stop("ppm should be positive")}
  if(rt < 0){stop("RT window should be positive")}
  if(!(toupper(RT.cor.method) %in% c("L", "P")) == TRUE) {stop("wrong RT correction method")}
  if(!(toupper(use.DB) %in% c("METABOLOMICS", "LIPIDOMICS")) == TRUE) {stop("wrong DB selected")}

  #(2) select DB
  if(toupper(use.DB) == "METABOLOMICS") {DB <- as.data.frame(sysdata$metabolite_db)}
  if(toupper(use.DB) == "LIPIDOMICS") {DB <- as.data.frame(sysdata$lipid_db)}


  if(!isTRUE(RT.cor) == TRUE) {
    if(dim(RT.DB)[1] != dim(DB)[1]) {stop("wrong RT correction DB, did you recalibrate your DB?")}
    if(toupper(RT.cor.method)  == "L") {DB$rt = RT.DB$rt_correct_l}
    if(toupper(RT.cor.method)  == "P") {DB$rt = RT.DB$rt_correct_p}
    }

  #(3) prepare the data, seperate MS1 and MS2
  pheno_levels <- levels(xset@phenoData$class)
  peak <- peakTable(xset)
  A = peak[, c(-1:-(7 + length(pheno_levels)))]
  B = t(A)
  class = row.names(xset@phenoData)
  mslevels = sub(".*(\\d+{1}).*$", "\\1", class)
  C <- cbind.data.frame(B, mslevels)
  peak_ms1 <- C[C$mslevels == "1",]
  peak_ms1$mslevels = NULL
  peak_ms1 <- t(peak_ms1)
  peak_ms2 <- C[C$mslevels == "2",]
  peak_ms2$mslevels = NULL
  peak_ms2 <- t(peak_ms2)
  peak_ms11 <- cbind(peak[, c(1:(7 + length(pheno_levels)))], peak_ms1)
  peak_ms22 <- cbind(peak[, c(1:(7 + length(pheno_levels)))], peak_ms2)
  mymz = peak_ms11$mz
  mymz = round(mymz, digits = 4)
  myrt = peak_ms11$rt
  myrt = round(myrt, digits = 2)

  #(4) search in database
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))
  Result <- vector("list", length(mymz))

  for (i in 1:length(mymz)) {
    ##(4.1) for MS1
    width <- options()$width * 0.3
    cat(paste0(rep(c(intToUtf8(0x2698), "="), i / length(mymz) * width), collapse = ''))
    cat(paste0(round(i / length(mymz) * 100), '% completed'))
    ## get [M+H]
    DB.list_H <- expand.grid.df(i, mymz[i], myrt[i], "[M+H]+", DB[, -c(6, 7, 8)])
    colnames(DB.list_H)[c(1:4, 9)] <- c("QueryID", "My.mz", "My.RT", "Adduct", "T.mz")
    ## get [M+Na]
    DB.list_Na <- expand.grid.df(i, mymz[i], myrt[i], "[M+Na]+", DB[, -c(5, 7, 8)])
    colnames(DB.list_Na)[c(1:4, 9)] <- c("QueryID", "My.mz", "My.RT", "Adduct", "T.mz")
    ## get [M+K]
    DB.list_K <- expand.grid.df(i, mymz[i], myrt[i], "[M+K]+", DB[, -c(5, 6, 8)])
    colnames(DB.list_K)[c(1:4, 9)] <- c("QueryID", "My.mz", "My.RT", "Adduct", "T.mz")

    ## get [M+NH4]
    DB.list_NH4 <- expand.grid.df(i, mymz[i], myrt[i], "[M+NH4]+", DB[, -c(5, 6, 7)])
    colnames(DB.list_NH4)[c(1:4, 9)] <- c("QueryID", "My.mz", "My.RT", "Adduct", "T.mz")
    ## combine them
    DB.list <- rbind(DB.list_H, DB.list_Na, DB.list_NH4)
    ## filter accoding to ppm and RT
    cal_ppm <- with(DB.list, (T.mz - My.mz) * 10^6 / T.mz)
    cal_ppm <-  round(cal_ppm, digits = 2)
    RT_diff <- with(DB.list, abs(My.RT - rt))
    RT_diff <- round(RT_diff, digits = 2)
    DB.list <- cbind(DB.list, Cal.ppm = cal_ppm, RT.dif = RT_diff)
    ## add a new column, My.msms
    ## rbind.data.frame does not work for the subsequent convertion of list tp DF, use rbind() here
    n.DB.list <- dim(DB.list)[1]
    DB.list <- cbind(DB.list, My.msms = rep(NA, n.DB.list))
    Result[[i]] = DB.list[(abs(DB.list$Cal.ppm) <= ppm) & DB.list$RT.dif <= rt, ]
    row.names(Result[[i]]) <- NULL

    ##(4.2) for MS2
    n_Result <- dim(Result[[i]])[1]
    if (n_Result > 0 & frag == TRUE){
      for (j in 1:n_Result) {
        ## prepare the data
        msms = as.numeric(strsplit(as.vector(Result[[i]][j, ]$msms.frgs), split = ";")[[1]])
        msms_rt = rep(Result[[i]][j, ]$rt, length(msms))
        T.frag <- cbind.data.frame(msms = msms, msms_rt = msms_rt)
        ms2 <- cbind.data.frame(ms2 = peak_ms22$mz, rt2 = peak_ms22$rt)
        ms2 <- round(ms2, digits = 4)
        myMS2 = expand.grid.df(ms2, T.frag)
        msms_ppm <- with(myMS2, (msms - ms2) * 10^6 / msms)
        msms_ppm <-  round(msms_ppm, digits = 2)
        msms_RT_diff <- with(myMS2, abs(rt2 - msms_rt))
        msms_RT_diff <- round(msms_RT_diff, digits = 2)
        myMS2 <- cbind(myMS2, Cal.ppm = msms_ppm, RT.dif = msms_RT_diff)
        getmsms <- myMS2[(abs(myMS2$Cal.ppm) <= ppm) & myMS2$RT.dif <= rt, ]
        ## check if msms can be found, and save them in the result
        if(dim(getmsms)[1] > 0) {
          Result[[i]][j, ]$My.msms <- paste(getmsms$ms2, collapse = ";")
        }
      }
    }

    if (i == length(mz)) cat(': Searching Done.')
    else cat('\014')
  }

  #(5) format the result
  search_result <- do.call(rbind.data.frame, Result)
  ## remove My.msms column when no MS/MS search is performed
  if(frag == FALSE) {search_result$My.msms = NULL}
  iden_id <- search_result$QueryID
  ## select identified rows in MS1
  iden_ms1 <- peak_ms11[iden_id, ]
  ## combine the search_result with identified MS1
  iden_result <- cbind.data.frame(search_result, iden_ms1[, -c(1:6)])
  non_iden_ms1 <- peak_ms11[-iden_id, ]
  iden_result$My.RT <- round(iden_result$My.RT/60, 2) # convert RT to min
  iden_result$rt <- round(iden_result$rt/60, 2)
  iden_result$RT.dif <- round(iden_result$RT.dif/60, 2)
  iden_result <- iden_result[order(iden_result$My.RT, decreasing = FALSE),] # order according to RT
  row.names(iden_result) <- NULL
  final_result <- list(Identified = iden_result, Not_Identified = non_iden_ms1)
  return(final_result)
}
