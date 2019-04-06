#' @title metabolite/lipid identification
#' @description metabolite/lipid identification using home database accoring to m/z and RT. MSMS fragments are also given.
#' @author Yonghui Dong
#' @param xset xcms object
#' @param RT.DB correction RT result from RT_correct()
#' @param use.DB which database is used for peak identification, choose between "Lipidomics" and "Metabolomics".
#' @param mode ionization mode, positive "+" or negative "-".
#' @param RT.cor should use corrected RT
#' @param RT.cor.method which corrected RT should be used? L: linear regression (default) or P: polynomial regression?
#' @param frag should search fragment
#' @param ppm mass tolerance, default value = 10
#' @param rt retention time search window, default value = 25
#' @param ms2.rm should remove MSMS data when it is included, default = TRUE
#' @param subgroup subset only specific groups for statistical test, default = NULL, which means statistics will be performed to all group pairs.
#' @param scale_group select groups needs to be scaled up or down, accoding to sample dilution or concentration.
#' @param scale_factor the scale factor, default = 1.
#' @importFrom xcms peakTable
#' @importFrom CAMERA annotate getPeaklist
#' @export
#'@examples
#'\dontrun{
#'my_compound <- what2(xset)
#'}

what2 <- function (xset, RT.DB = NULL, use.DB = "METABOLOMICS", mode = "+", RT.cor = TRUE,
                   RT.cor.method = "L", frag = T, ppm = 10, rt = 25, ms2.rm = T, subgroup = NULL,
                   scale_group = NULL, scale_factor = 1) {

  cat("\n(1) Checking input parameters...");
  #(1) input check
  pheno_levels <- levels(xset@phenoData$class)
  if(class(xset) != "xcmsSet") {stop("the input object is not an xcmsSet object")}
  if(!is.numeric(ppm)){stop("invalid calss of ppm threshold: not numeric")}
  if(!is.numeric(rt)){stop("invalid calss of RT threshold: not numeric")}
  if(ppm < 0){stop("ppm should be positive")}
  if(rt < 0){stop("RT window should be positive")}
  if(!(toupper(RT.cor.method) %in% c("L", "P")) == TRUE) {stop("wrong RT correction method")}
  if(!(toupper(use.DB) %in% c("METABOLOMICS", "LIPIDOMICS")) == TRUE) {stop("wrong DB selected")}
  if(!(mode %in% c("+", "-")) == TRUE) {stop("wrong ion mode, select '+' or '-'")}
  if(is.null(subgroup) == FALSE & all(subgroup %in% pheno_levels) == FALSE)
  {stop("selected subgroup(s) do not exist in your data")}
  if(is.null(scale_group) == FALSE & all(scale_group %in% pheno_levels) == FALSE)
  {stop("selected scale group(s) do not exist in your data")}
  if(!is.numeric(scale_factor)){stop("invalid calss of scale factor: not numeric")}
  if(scale_factor < 0){stop("scale factor should be positive")}

  #(2) select DB
  if(toupper(use.DB) == "METABOLOMICS" & mode == "+") {DB <- as.data.frame(sysdata$Metabolite_DB_pos)}
  if(toupper(use.DB) == "METABOLOMICS" & mode == "-") {DB <- as.data.frame(sysdata$Metabolite_DB_neg)}
  if(toupper(use.DB) == "LIPIDOMICS" & mode == "+") {DB <- as.data.frame(sysdata$Lipid_DB_pos)}

  if(isTRUE(RT.cor) == TRUE) {
    if(dim(RT.DB)[1] != dim(DB)[1]) {stop("wrong RT correction DB, did you recalibrate your DB?")}
    if(toupper(RT.cor.method)  == "L") {DB$T.RT = RT.DB$rt_correct_l}
    if(toupper(RT.cor.method)  == "P") {DB$T.RT = RT.DB$rt_correct_p}
  }

  cat("passed");

  #(3) prepare the data, seperate MS1 and MS2
  cat("\n(2) Deisotoping...");
  pheno_levels <- levels(xset@phenoData$class)
  peak <- peakTable(xset)

  ##(3.1) add scaling facors. for instance, some samples were diluted 10 times, so the intensity should be multiplied by 10
  my_meta <- xset@phenoData
  if (is.null(scale_group) == FALSE){
    get_cnames <- row.names(my_meta)[my_meta$class %in% scale_group]
    peak[, (colnames(peak) %in% get_cnames)] <- peak[, (colnames(peak) %in% get_cnames)] * scale_factor
  }

  ##(3.2) change colnames to ease futher data analysis
  colnames(peak)[-c(1:(7 + length(pheno_levels)))] <- paste(my_meta$class, "_", row.names(my_meta), sep = "")

  ##(3.3) deisotoping
  anI <- annotate(xset, ppm = 10, multiplier = 2, quick = TRUE, calcIso = TRUE)
  iso_peaklist <- getPeaklist(anI)
  iso_peaklist$isotopes <- sub("\\[.*?\\]", "", iso_peaklist$isotopes)
  peak <- peak[iso_peaklist$isotopes == '' | iso_peaklist$isotopes == '[M]+', ]

  ##(3.4) prepare the data
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
  cat("Done");

  #(4) get fold change and p-values

  cat("\n(3) Calculating fold change and performing pair-wise statistical test...");
  data_summary <- mysummary(xset, ms2.rm = ms2.rm, subgroup = subgroup,
                            scale_group = scale_group, scale_factor = scale_factor)
  ## deisotoping
  data_summary <- data_summary[iso_peaklist$isotopes == '' | iso_peaklist$isotopes == '[M]+', ]

  ## get the sample name which has the max intensity of each m/z. It is useful for further manual identification
  ## only from ms1
  tmp1 <- peak_ms11[, c(-1:-(7 + length(pheno_levels)))]
  sample_index <- as.matrix(apply(tmp1, 1, which.max))
  sample_index <- as.vector(sample_index, mode = "numeric")
  sample_name <- colnames(tmp1)[sample_index]
  data_summary <- cbind(Max_Sample = sample_name, tmp1, data_summary)
  cat("done");

  #(5) search in database
  cat("\n(4) Searching against database...");
  ## supress warnings
  options(warn=-1)

  expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))
  Result <- vector("list", length(mymz))

  for (i in 1:length(mymz)) {
    ##(4.1) for MS1
    width <- options()$width * 0.3
    cat(paste0(rep(c(intToUtf8(0x2698), "="), i / length(mymz) * width), collapse = ''))
    cat(paste0(round(i / length(mymz) * 100), '% Database searching completed...'))
    DB.list <- expand.grid.df(i, mymz[i], myrt[i], DB)
    colnames(DB.list)[c(1:3)] <- c("QueryID", "My.mz", "My.RT")

    ## filter accoding to ppm and RT
    cal_ppm <- with(DB.list, (T.mz - My.mz) * 10^6 / T.mz)
    cal_ppm <-  round(cal_ppm, digits = 2)
    RT_diff <- with(DB.list, abs(My.RT - T.RT))
    RT_diff <- round(RT_diff, digits = 2)
    DB.list <- cbind(DB.list, Cal.ppm = cal_ppm, RT.dif = RT_diff)

    ## add a new column, My.fragment
    ## rbind.data.frame does not work for the subsequent convertion of list tp DF, use rbind() here
    n.DB.list <- dim(DB.list)[1]
    DB.list <- cbind(DB.list, My.msms = rep(NA, n.DB.list))
    Result[[i]] = DB.list[(abs(DB.list$Cal.ppm) <= ppm) & DB.list$RT.dif <= rt, ]
    row.names(Result[[i]]) <- NULL

    ##(4.2) for MS2
    n_Result <- dim(Result[[i]])[1]
    if (n_Result > 0 & frag == TRUE){
      for (j in 1:n_Result) {
        ## check if msms exist, if no, skip this iteration
        if(all(is.na(Result[[i]][j, ]$Fragment)) == TRUE) next
        ## prepare the data
        msms = as.numeric(strsplit(as.vector(Result[[i]][j, ]$Fragment), split = ";")[[1]])
        msms_rt = rep(Result[[i]][j, ]$T.RT, length(msms))
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
  cat("\n(5) Formating the result...");
  search_result <- do.call(rbind.data.frame, Result)
  search_result$rt = NULL
  ## remove My.msms column when no MS/MS search is performed
  if(frag == FALSE) {search_result$My.msms = NULL}
  iden_id <- search_result$QueryID
  ## select identified rows in MS1
  iden_ms1 <- data_summary[iden_id, ]

  ## combine the search_result with identified MS1
  iden_result <- cbind.data.frame(search_result, iden_ms1)
  non_iden_ms1 <- data_summary[-iden_id, ]
  # add mz, RT infor
  non_iden_ms1 <- cbind.data.frame(peak[-iden_id, c(1:(7 + length(pheno_levels)))], non_iden_ms1)
  iden_result$My.RT <- round(iden_result$My.RT/60, 2) # convert RT to min
  iden_result$T.RT <- round(iden_result$T.RT/60, 2)
  iden_result$RT.dif <- round(iden_result$RT.dif/60, 2)
  iden_result <- iden_result[order(iden_result$My.RT, decreasing = FALSE),] # order according to RT
  row.names(iden_result) <- NULL
  final_result <- list(Full_peak = peak, Identified = iden_result, Not_Identified = non_iden_ms1)
  cat("done");
  return(final_result)
}
