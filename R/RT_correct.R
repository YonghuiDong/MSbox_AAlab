#' @title Retention time calibration
#' @description calibrate RT according to reference peaks
#' @author Yonghui Dong
#' @param xset xcms object
#' @param ppm mass tolerance, default value = 10
#' @param rt retention time search window, default value = 50
#' @param anchors lipidomics or metabolomics reference peaks
#' @importFrom stats lm prcomp predict sd
#' @export
#' @examples
#'  a = 1+1



RT_correction <- function(xset, ppm = 10, rt = 50, anchors = "LIPIDOMICS"){

  #(2) define anchors and select DB
  pre_anchors = c("LysoPC 18:1(1)", "lysoPC 16:0", "lysoPC 16:0", "lysoPC 18:2(1)", "PC 34:3","PC 36:2","PC 36:3 (1)",
              "PC 36:4 (1)","PC 36:5","PC 36:6","LysoPC 18:0","lysoPC 18:3(1)","PE 34:2","PE 36:5","TAG 54:5",
              "GlcCer t18:1/h24:0", "TAG 52:0", "TAG 58:3")
  if(anchors == "LIPIDOMICS") {anchors <- pre_anchors}
  if(anchors == "LIPIDOMICS") {DB <- as.data.frame(sysdata$lipid_db)}

  #(3) pre-annotation for RT correction

  ##(3.1) prepare the data, seperate MS1 and MS2
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

  ##(3.2) search in database
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))
  Result <- vector("list", length(mymz))

  for (i in 1:length(mymz)) {
    ## for MS1
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
    Result[[i]] = DB.list[(abs(DB.list$Cal.ppm) <= ppm) & DB.list$RT.dif <= rt, ]
    row.names(Result[[i]]) <- NULL
    if (i == length(mz)) cat(': Searching Done.')
    else cat('\014')
  }

  #(4) format the result
  my_iden <- do.call(rbind.data.frame, Result)


  ##(5) calibration
  my.anchors = my_iden[my_iden$compound %in% anchors,]$My.RT
  if(length(my.anchors) == 0) {stop("No anchors were found!")}
  anchors.library = my_iden[my_iden$compound %in% anchors,]$rt
  unique.anchors = length(unique(my_iden[my_iden$compound %in% anchors,]$compound))
  RT.range <- round(range(my.anchors)/60, 2)
  cat('\014')
  cat(paste(unique.anchors, "reference anchors are used for RT correction. The RT range is between", RT.range[1], "-", RT.range[2], "min"))

  ## linear may be safer than polynomial (overfitting problem)

  model = lm(my.anchors ~ anchors.library)
  modelp <- lm(my.anchors ~ poly(anchors.library, 3))

  adj.R = round(summary(model)$adj.r.squared, 3)
  adj.R_p = round(summary(modelp)$adj.r.squared, 3)

  predicted.rt <- as.data.frame(predict(model, newdata = data.frame(anchors.library = DB$rt),
                                        interval='confidence', level=0.99))
  predicted.rt_p <- as.data.frame(predict(modelp, newdata = data.frame(anchors.library = DB$rt),
                                          interval='confidence', level=0.99))

  par(mfrow=c(2,1))
  plot(x = DB$rt,  y = predicted.rt$fit, col = "red", pch = 16, xlab = "Original RT", ylab = "Predicted RT")
  abline(lm(my.anchors ~ anchors.library), col = "blue", lwd=3)
  legend("topleft", paste("adj.R^2 = ", adj.R, sep =""), bty = "n")
  title("Linear Regression")

  plot(x = DB$rt,  y = predicted.rt_p$fit, col = "red", pch = 16, xlab = "Original RT", ylab = "Predicted RT")
  abline(lm(my.anchors ~ anchors.library), col = "blue", lwd=3)
  legend("topleft", paste("adj.R^2 = ", adj.R_p, sep =""), bty = "n")
  title("Polynomia Regression, degree = 3")

  ## replace DB RT with new fitted RT
  DB$rt_correct_l = predicted.rt$fit
  DB$rt_correct_p = predicted.rt_p$fit

  cat(sep="\n\n")
  cat(sep="\n\n")
  cat("Attention:\n")
  cat(sep="\n\n")
  cat("(1) For perfect RT correction, you would have Predicted RT = Original RT.\n")
  cat("(2) If there is a strong distortion, please carefully re-select your anchors, and re-do RT correction. \n")
  cat("(3) Linear and polynomial regresstion are used for RT correction. Select the prefered method for metabolite identification. The default one is linear regression. \n")
  return(DB)
  }




