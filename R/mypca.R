#' @title fast PCA
#' @description perform PCA from xcms object
#' @author Yonghui Dong
#' @param xset xcms object
#' @param centering centering, default = TRUE
#' @param scaling scaling method, default is scaling = "none". You can choose "auto" or "pareto"
#' @param ms2.rm remove MSMS data when it is included.
#' @param x PCA X axis, default is PC1
#' @param y PCA Y axis, defult is PC2
#' @param exclude exclude some classes of samples
#' @param size dot size
#' @param interactive should interactive figure be plotted? default = TRUE. If you want to save the result
#' in high resolution, use non interative plot.
#' @param scale_group select groups needs to be scaled.
#' @param scale_factor the scale factor, default = 1.
#' @param ... other parameters
#' @importFrom ggplot2 autoplot theme_bw
#' @import ggfortify
#' @importFrom graphics abline legend par plot title
#' @importFrom xcms peakTable
#' @importFrom plotly ggplotly
#' @export
#'@examples
#'\dontrun{
#'getpca <- mypca(xset)
#'}

mypca <- function(xset, centering = T, scaling = "none", ms2.rm = FALSE, x = 1, y = 2,
                 size = 1.5, exclude = NULL, scale_group = NULL, scale_factor = 1,
                 interactive = T, ...) {
  #(1) check input
  ## check object type
  if(class(xset) != "xcmsSet") {stop("the input object is not an xcmsSet object")}
  ## check xset phenoData
  pheno_levels <- levels(xset@phenoData$class)
  if(length(pheno_levels) < 2) {stop("at least two sample groups should be included")}
  ## check scaling
  if(!(toupper(scaling) %in% c("PARETO", "AUTO", "NONE")) == TRUE) {stop("wrong scaling method")}

  #(2) extract xcms information
  peak <- peakTable(xset)

  ##(2.1) add scaling facors. for instance, some samples were diluted 10 times, so the intensity should be multiplied by 10
  my_meta <- xset@phenoData
  if (is.null(scale_group) == FALSE){
    get_cnames <- row.names(my_meta)[my_meta$class %in% scale_group]
    peak[, (colnames(peak) %in% get_cnames)] <- peak[, (colnames(peak) %in% get_cnames)] * scale_factor
  }

  ##(2.1) add rownames
  mz = round(peak$mz, 3)
  rt = round(peak$rt/60, 2)
  mynames <- paste(mz, "_", rt, sep ="")
  row.names(peak) <- mynames
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

  Group <- B$Group
  B1 <- B[, -which(names(B) %in% "Group")]
  B1 = B1[,apply(B1,2,function(x) !all(x==0))]

  #(3) apply scaling
  Group <- B$Group
  if(toupper(scaling) == "PARETO") {
    B2 <- apply(B1, 2, function(x) x/sqrt(sd(x)))
  } else if(toupper(scaling) == "AUTO") {
    B2 <- scale(B1, center = FALSE, scale = TRUE)
  } else if(toupper(scaling) == "NONE") {
    B2 = B1
  }

  #(4) excluding
  ## only select subgroups if subgroup is not NULL
  D = cbind.data.frame(B2, Group)
  if(is.null(exclude) == TRUE) {
    D_new = D
  } else {
    D_new <- D[!(toupper(D$Group) %in% toupper(exclude)),]
    ## drop factors
    D_new$Group <- factor(D_new$Group)
  }

  #(5) perform PCA
  B3 = D_new
  B3$Group = NULL

  p <- autoplot(prcomp(B3, center = centering), D_new, colour = "Group",
                size = size, x = x, y = y, ...) +
    theme_bw()
  p2 <- ggplotly(p)

  ## choose to use interactive or non-interactive plot
  if(isTRUE(interactive) == TRUE) {
    return(p2)
  } else {
    return(p)
  }

}
