#' @title Prefilter
#' @description prefiltering isotopically labeled analytes according to the experiment design.
#' @param xset xcms object.
#' @param ms2.rm remove MSMS data when it is included, default = TRUE
#' @param subgroup subset the xcms groups. The name should be the same as in phboData$class. default = NULL, which means no subset will be performed.
#' @importFrom xcms peakTable
#' @importFrom stats lm
#' @importFrom stats aov
#' @importFrom stats TukeyHSD
#' @importFrom xcms peakTable
#' @return a filtered peaklist
#' @export
#' @examples
#'\dontrun{
#' my_summary <- mysummary(xset)
#'}

mysummary <- function(xset, ms2.rm = ms2.rm, subgroup = subgroup){

  #(1) check input
  ##(1.1) check xset phenoData
  pheno_levels <- levels(xset@phenoData$class)
  ##(1.2) check subsetgroup
  if(is.null(subgroup) == FALSE & all(subgroup %in% pheno_levels) == FALSE)
  {stop("selected subgroup(s) do not exist in your data")}

  #(2) extract xcms information
  peak <- peakTable(xset)
  ##(2.1) deisotoping
  an <- xsAnnotate(xset)
  anG <- groupFWHM(an)
  anI <- findIsotopes(an)
  iso_peaklist <- getPeaklist(anI)
  iso_peaklist$isotopes <- sub("\\[.*?\\]", "", iso_peaklist$isotopes)
  peak <- peak[iso_peaklist$isotopes == '' | iso_peaklist$isotopes == '[M]+', ]
  ##(2.2) prepare the data
  A = peak[, c(-1:-(7 + length(pheno_levels)))]
  B = cbind.data.frame(t(A), Group = xset@phenoData$class)

  #(3) remove MS2 if exist
  if (isTRUE(ms2.rm) == TRUE) {
    class = row.names(xset@phenoData)
    mslevels = sub(".*(\\d+{1}).*$", "\\1", class)
    C <- cbind.data.frame(B, mslevels)
    peak_ms1 <- C[C$mslevels == "1",]
    peak_ms1$mslevels = NULL
    B <- peak_ms1
  }

  #(4) only select subgroups if subgroup is not NULL
  if(is.null(subgroup) == TRUE) {
    peaklist_new = B
  } else {
    peaklist_new <- B[(B$Group %in% subgroup),]
    ## drop factors
    peaklist_new$Group <- factor(peaklist_new$Group)
  }

  ##(5) calculating fold change
  fold <- fold(peaklist_new)

  ##(6) statistical test
  stat <- getp(peaklist_new)

  ##(7) combind result
  my_summary <- cbind.data.frame(A, fold, stat)
  return(my_summary)
}

