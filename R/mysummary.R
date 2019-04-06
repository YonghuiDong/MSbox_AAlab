#' @title Prefilter
#' @description prefiltering isotopically labeled analytes according to the experiment design.
#' @param xset xcms object.
#' @param ms2.rm should remove MSMS data when it is included, default = TRUE
#' @param subgroup subset only specific groups for statistical test, default = NULL, which means statistics will be performed to all group pairs.
#' @param scale_group select groups needs to be scaled up or down, accoding to sample dilution or concentration.
#' @param scale_factor the scale factor, default = 1.
#' @importFrom xcms peakTable
#' @importFrom stats lm
#' @importFrom stats aov
#' @importFrom stats TukeyHSD
#' @importFrom xcms peakTable
#' @importFrom CAMERA annotate getPeaklist
#' @return a filtered peaklist
#' @export
#' @examples
#'\dontrun{
#' my_summary <- mysummary(xset)
#'}

mysummary <- function(xset, ms2.rm, subgroup, scale_group, scale_factor){

  #(1) extract xcms information
  peak <- peakTable(xset)
  pheno_levels <- levels(xset@phenoData$class)
  ##(1.1) add scaling facors. e.g., some groups were diluted 10 times before analysis, so the intensity should be multiplied by 10
  my_meta <- xset@phenoData
  if (is.null(scale_group) == FALSE){
    get_cnames <- row.names(my_meta)[my_meta$class %in% scale_group]
    peak[, (colnames(peak) %in% get_cnames)] <- peak[, (colnames(peak) %in% get_cnames)] * scale_factor
  }
  ##(1.2) change colnames to ease futher data analysis
  colnames(peak)[-c(1:(7 + length(pheno_levels)))] <- paste(my_meta$class, "_", row.names(my_meta), sep = "")

  ##(1.3) prepare the data
  A = peak[, c(-1:-(7 + length(pheno_levels)))]
  B = cbind.data.frame(t(A), Group = xset@phenoData$class)

  #(2) remove MS2 if exist
  if (isTRUE(ms2.rm) == TRUE) {
    class = row.names(xset@phenoData)
    mslevels = sub(".*(\\d+{1}).*$", "\\1", class)
    C <- cbind.data.frame(B, mslevels)
    peak_ms1 <- C[C$mslevels == "1",]
    peak_ms1$mslevels = NULL
    B <- peak_ms1
  }

  #(3) only select subgroups if subgroup is not NULL
  if(is.null(subgroup) == TRUE) {
    peaklist_new = B
  } else {
    peaklist_new <- B[(B$Group %in% subgroup),]
    ## drop factors
    peaklist_new$Group <- factor(peaklist_new$Group)
  }

  ##(4) calculating fold change
  fold <- fold(peaklist_new)

  ##(5) statistical test
  stat <- getp(peaklist_new)

  ##(6) combind result
  my_summary <- cbind.data.frame(fold, stat)
  return(my_summary)
}

