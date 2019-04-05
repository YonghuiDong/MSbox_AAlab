#' @title Retention time calibration
#' @description calibrate RT according to reference peaks
#' @author Yonghui Dong
#' @param xset xcms object
#' @param ppm mass tolerance, default value = 10
#' @param rt retention time search window (in second), default value = 50 s.
#' @param anchors lipidomics or metabolomics reference peaks
#' @param use.DB which database is used for peak identification, choose between "Lipidomics" and "Metabolomics".
#' @param mode ionization mode, positive "+" or negative "-".
#' @importFrom stats lm prcomp predict sd median
#' @importFrom xcms peakTable
#' @importFrom CAMERA annotate getPeaklist
#' @export
#' @examples
#'\dontrun{
#'RT.DB <- RT_correction(xset, use.DB = "Lipidomics", mode = "+")
#'}



RT_correction <- function(xset, ppm = 10, rt = 50, anchors = NULL, use.DB = "LIPIDOMICS", mode = "+"){

  #(1) check input
  if(class(xset) != "xcmsSet") {stop("the input object is not an xcmsSet object")}
  if(!is.numeric(ppm)){stop("invalid calss of ppm threshold: not numeric")}
  if(!is.numeric(rt)){stop("invalid calss of RT threshold: not numeric")}
  if(ppm < 0){stop("ppm should be positive")}
  if(rt < 0){stop("RT window should be positive")}
  if(!(mode %in% c("+", "-")) == TRUE) {stop("wrong ion mode, select '+' or '-'")}


  #(2) define anchors and select DB
  ##(2,1) for lipidomics pos
  pre_anchors_lipid_pos = c("LysoPC 18:1(1)", "lysoPC 16:0", "lysoPC 16:0", "lysoPC 18:2(1)", "PC 34:3","PC 36:2","PC 36:3 (1)",
              "PC 36:4 (1)","PC 36:5","PC 36:6","LysoPC 18:0","lysoPC 18:3(1)","PE 34:2","PE 36:5","TAG 54:5",
              "GlcCer t18:1/h24:0", "TAG 52:0", "TAG 58:3")
  if(toupper(use.DB) == "LIPIDOMICS" & mode == "+" & is.null(anchors) == TRUE) {anchors <- pre_anchors_lipid_pos}
  if(toupper(use.DB) == "LIPIDOMICS" & mode == "+") {DB <- as.data.frame(sysdata$Lipid_DB_pos)}

  ##(2.2) for lipidomics neg (no DB now)

  ##(2.3) for metabolomics pos (anchors need to be optimized, not ready)
  pre_anchors_metabolite_pos = c("Tyrosine (S)", "Phenylalanine (S)", "Caffeoylputrescine", "Tryptophan (S)",
                                 "Coumaric acid hexose", "Quercetin-dihexose", "p-Coumaric acid", "Ferulic acid",
                                 "Sinapic acid", "Kaempferol hexose", "Dehydrotomatine II", "Dehydrotomatine III",
                                 "a-Tomatine (S)", "Tomatidine", "Naringenin hexose", "Rutin (S)", "Tomatidine isomer",
                                 "Tomato-UGA 3", "Coumaroyl diferuoyl spremidine")
  if(toupper(use.DB) == "METABOLOMICS" & mode == "+" & is.null(anchors) == TRUE) {anchors <- pre_anchors_metabolite_pos}
  if(toupper(use.DB) == "METABOLOMICS" & mode == "+") {DB <- as.data.frame(sysdata$Metabolite_DB_pos)}

  ##(2.4) for metabolomics neg (anchors need to be optimized, not ready)
  pre_anchors_metabolite_neg = c("Malic acid Isomer 1", "Tyrosine (S)","p-Coumaric acid (S)", "Galloyl quinic acid",
                                 "Phenylalanine (S)", "Di-hydroxybenzoic acid-Hexose", "Phenylalanine derivative",
                                 "Neochlorogenic acid", "Caffeic acid-hexose Isomer 1", "Caffeic acid", "Tryptophan (S)",
                                 "Kaempferol 3-O-triglucoside-7-O-glucoside", "procyanidin_dimer", "Chlorogenic acid",
                                 "Quercetin-dihexose-pentose-deoxyhexose", "Quercetin-3-O-(caffeoyl) trihexose",
                                 "1-O-Caffeoylquinic acid", "Coumaric acid-Hexose Isomer 3", "Ferulic acid Isomer 1",
                                 "Quercetin hexose-hexose", "Quercetin-dihexose-deoxyhexose", "Rutin",
                                 "Dehydrotomatine (S)", "Tomatine (S)", "Labdanolic acid", "Grandiflorenic acid", "Stearolic acid",
                                 "Colladonin", "Hocogenin", "Malic acid Isomer 2", "p-Coumaric acid (S)", "Feruloyl quinic acid",
                                 "Quercetin hexose-hexose", "Quercetin-Hex-deoxyHex", "quercetin 3-caffeoyl dihexose", "Kaempferol glucuronide hexose",
                                 "Ferulic acid", "Kaempferol dihexose", "Kaempferol-hexose-deoxyhexose Isomer 1",
                                 "Sinapic acid (S)", "quercetin-3-O-feruloyl dihexose Isomer 2", "Esculeoside B",
                                 "Coumaroyl malate Isomer 1", "kaempferol-3-O-caffeoyl sophorotrioside", "Quercetin glucuronide",
                                 "Malic acid Isomer 3", "Kaempferol hexose", "Quercetin deoxyhexose", "isorhamnetin hexose",
                                 "Naringenin hexose", "Isorhamnetin hexose", "Kaempferol deoxyhexose", "Quercetin deoxyhexose",
                                 "Tomato-UGA 15", "feruloylsinapoylglucose", "Naringenin (S)")

  if(toupper(use.DB) == "METABOLOMICS" & mode == "-" & is.null(anchors) == TRUE) {anchors <- pre_anchors_metabolite_neg}
  if(toupper(use.DB) == "METABOLOMICS" & mode == "-") {DB <- as.data.frame(sysdata$Metabolite_DB_neg)}


  #(3) pre-annotation for RT correction

  ##(3.1) prepare the data, seperate MS1 and MS2
  pheno_levels <- levels(xset@phenoData$class)
  peak <- peakTable(xset)

  ##(3.2) deisotoping
  anI <- annotate(xset, ppm = 10, quick = TRUE, calcIso = TRUE)
  iso_peaklist <- getPeaklist(anI)
  iso_peaklist$isotopes <- sub("\\[.*?\\]", "", iso_peaklist$isotopes)
  peak <- peak[iso_peaklist$isotopes == '' | iso_peaklist$isotopes == '[M]+', ]

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

    ##(4.1) for MS1
    width <- options()$width * 0.3
    cat(paste0(rep(c(intToUtf8(0x2698), "="), i / length(mymz) * width), collapse = ''))
    cat(paste0(round(i / length(mymz) * 100), '% completed'))
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

    if (i == length(mz)) cat(': Searching Done.')
    else cat('\014')
  }

  #(4) format the result
  my_iden <- do.call(rbind.data.frame, Result)

  ##(5) calibration
  my.anchors <- my_iden[my_iden$Compound %in% anchors,]$My.RT
  my.compound <- my_iden[my_iden$Compound %in% anchors,]$Compound
  if(length(my.anchors) == 0) {stop("No anchors were found!")}

  ## sometimes more than 1 compound were annotated as the same compound, then their median rt are used for RT correction
  my.anchors2 <- sapply(unique(my.compound), function(x) median(my.anchors[my.compound == x]))
  anchors.library <- unique(my_iden[my_iden$Compound %in% anchors,]$T.RT)


  ## provide anchor information
  anchor.data <- cbind.data.frame(compound = unique(my.compound), RT = round(my.anchors2/60, 2))
  anchor.data <- anchor.data[order(anchor.data$RT, decreasing = F), ]
  unique.anchors = dim(anchor.data)[1]
  RT.range <- round(range(anchor.data$RT), 2)
  cat('\014')
  cat(paste(unique.anchors, "references are used for RT correction. RT range is between", RT.range[1], "-", RT.range[2], "min"))
  cat(sep="\n\n")
  cat(paste0(rep("=",  width), collapse = ''))
  cat(sep="\n\n")
  print(anchor.data)
  cat(paste0(rep("=",  width), collapse = ''))


  ## linear may be safer than polynomial (overfitting problem)
  model = lm(my.anchors2 ~ anchors.library)
  modelp <- lm(my.anchors2 ~ poly(anchors.library, 3))

  adj.R = round(summary(model)$adj.r.squared, 3)
  adj.R_p = round(summary(modelp)$adj.r.squared, 3)

  predicted.rt <- as.data.frame(predict(model, newdata = data.frame(anchors.library = DB$T.RT),
                                        interval='confidence', level=0.99))
  predicted.rt_p <- as.data.frame(predict(modelp, newdata = data.frame(anchors.library = DB$T.RT),
                                          interval='confidence', level=0.99))

  par(mfrow=c(2,1))
  plot(x = DB$T.RT,  y = predicted.rt$fit, col = "red", pch = 16, xlab = "Original RT (s)", ylab = "Predicted RT (s)")
  abline(lm(my.anchors2 ~ anchors.library), col = "blue", lwd=3)
  legend("topleft", paste("adj.R^2 = ", adj.R, sep =""), bty = "n")
  title("Linear Regression")

  plot(x = DB$T.RT,  y = predicted.rt_p$fit, col = "red", pch = 16, xlab = "Original RT (s)", ylab = "Predicted RT (s)")
  abline(lm(my.anchors2 ~ anchors.library), col = "blue", lwd=3)
  legend("topleft", paste("adj.R^2 = ", adj.R_p, sep =""), bty = "n")
  title("Polynomia Regression, degree = 3")

  ## replace DB RT with new fitted RT
  DB$rt_correct_l = predicted.rt$fit
  DB$rt_correct_p = predicted.rt_p$fit

  cat(sep="\n\n")
  cat(sep="\n\n")
  cat("Attention:\n")
  cat(sep="\n\n")
  cat("Linear and polynomial regresstion are used for RT correction. Select the prefered method for metabolite identification. The default one is linear regression. \n")
  return(DB)
  }




