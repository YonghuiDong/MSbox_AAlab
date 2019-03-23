#' @title volcano plot
#' @description volcano plot.
#' @author Yonghui Dong
#' @param result result from what2() function
#' @param compare_group which groups you want to compare, i.e. c("WT", "Treat1"), only two groups are allowed
#' @param FC select fold change values, default = 2
#' @param p select p value, default = 0.05
#' @param show_plot volcano plot should be plotted according to "iden", identified compounds or "non_iden", non identified compounds. default = "iden"
#' @param interactive should interactive plot be shown? default = TRUE
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @export
#'@examples
#'\dontrun{
#'myvolcano(result, compare_group = c("WT", "JA"))
#'}

myvolcano <- function(result, compare_group, FC = 2, p = 0.05, show_plot = "iden", interactive = TRUE) {

  #(1)check input
  if(!is.numeric(FC)){stop("invalid calss of fold change: not numeric")}
  if(!is.numeric(p)){stop("invalid calss of p value: not numeric")}
  if(p > 1){stop("p value should between 0 and 1")}
  ## supress "no visible binding for global variable" notes in ggplot2
  Compoud_QueryID <- Fold_change <- aes <- my_cor <- p.value <- NULL

  #(2) get index
  F_id <- paste("Fold_", compare_group[1], "_VS_", compare_group[2], sep = "")
  P_id1 <- paste("P_", compare_group[1], "_VS_", compare_group[2], sep = "")
  P_id2 <- paste("P_", compare_group[2], "_VS_", compare_group[1], sep = "")

  #(3) plot

  if(toupper(show_plot) == "IDEN") {
    ##3.1 for identified result
    iden <- result$Identified
    ## select which p value colname exist.
    if ((P_id1 %in% colnames(iden)) == TRUE) {
      P_id = P_id1
    } else if ((P_id2 %in% colnames(iden)) == TRUE) {
      P_id = P_id2
    } else {stop("The selected groups do not exist")}
    ## prepare data frame for plot
    F_iden <- iden[, F_id]
    P_iden <- iden[, P_id]
    new_iden <- cbind.data.frame(Compoud_QueryID = iden$QueryID,
                                 Fold_change = log2(F_iden),
                                 p.value = -log10(P_iden))

    ## sometimes there are inf, NAN values in FC and p-value, so i replace these values with max(FC) and 0 for FC and p, respectively
    if(length(new_iden[!is.finite(new_iden$Fold_change), ]$Fold_change) > 0)
    {new_iden[!is.finite(new_iden$Fold_change), ]$Fold_change <- max(abs(new_iden[is.finite(new_iden$Fold_change),]$Fold_change))}
    if(length(new_iden[!is.finite(new_iden$p.value), ]$p.value) > 0)
    {new_iden[!is.finite(new_iden$p.value), ]$p.value <- 0}

    ## add color information
    new_iden$my_cor = "NS"
    con1 <- new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value < -log10(p)),]
    if(dim(con1)[1] > 0)
    {new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value < -log10(p)),]$my_cor =
      paste("|FC|", ">=", FC, "&", "p", ">", p, sep = " ")
    }

    con2 <- new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value >= -log10(p)),]
    if(dim(con2)[1] > 0)
    {
      new_iden[(abs(new_iden$Fold_change) >= log2(FC) & new_iden$p.value >= -log10(p)),]$my_cor =
        paste("|FC|", ">=", FC, "&", "p", "<=", p, sep = " ")
    }

    con3 <- new_iden[(abs(new_iden$Fold_change) < log2(FC) & new_iden$p.value >= -log10(p)),]
    if(dim(con3)[1] > 0)
    {
      new_iden[(abs(new_iden$Fold_change) < log2(FC) & new_iden$p.value >= -log10(p)),]$my_cor =
        paste("|FC|", "<", FC, "&", "p", "<=", p, sep = " ")
    }

    ##plot
    myplot <- ggplot(data = new_iden, aes(x = Fold_change, y = p.value, label = Compoud_QueryID, col = my_cor)) +
      geom_point(alpha = 1, size = 1) +
      geom_text(hjust=0, vjust=0, alpha = 1) +
      scale_color_manual(values=c("steelblue1", "Tomato", "springgreen3", "Grey")) +
      theme_bw() +
      geom_vline(xintercept = log2(FC), col = "black", linetype = "dotted", size = 0.6) +
      geom_vline(xintercept = -log2(FC), col = "black", linetype = "dotted", size = 0.7) +
      geom_hline(yintercept = -log10(p), col = "black", linetype = "dotted", size = 0.6) +
      ggtitle("Identified compound") +
      xlab("Log(2) fold change") +
      ylab("-Log(10) p-value") +
      labs(color = "")
  }

  if(toupper(show_plot) =="NON_IDEN") {
    ##3.2 for non-identified result
    non_iden <- result$Not_Identified
    ## select which p value colname exist.
    if ((P_id1 %in% colnames(non_iden)) == TRUE) {
      P_id = P_id1
    } else if ((P_id2 %in% colnames(non_iden)) == TRUE) {
      P_id = P_id2
    } else {stop("The selected groups do not exist")}
    ## prepare the data frame
    F_niden <- non_iden[, F_id]
    P_niden <- non_iden[, P_id]
    new_niden <- cbind.data.frame(Compoud_QueryID = row.names(non_iden),
                                  Fold_change = log2(F_niden),
                                  p.value = -log10(P_niden))
    ## sometimes there are inf, NAN values in FC and p-value, so i replace these values with max(FC) and 0 for FC and p, respectively
    if(length(new_niden[!is.finite(new_niden$Fold_change), ]$Fold_change) > 0)
    {new_niden[!is.finite(new_niden$Fold_change), ]$Fold_change <- max(abs(new_niden[is.finite(new_niden$Fold_change),]$Fold_change))}
    if(length(new_niden[!is.finite(new_niden$p.value), ]$p.value) > 0)
    {new_niden[!is.finite(new_niden$p.value), ]$p.value <- 0}
    ## add color information
    new_niden$my_cor = "NS"
    con1 <- new_niden[(abs(new_niden$Fold_change) >= log2(FC) & new_niden$p.value < -log10(p)),]
    if(dim(con1)[1] > 0){
      new_niden[(abs(new_niden$Fold_change) >= log2(FC) & new_niden$p.value < -log10(p)),]$my_cor =
        paste("|FC|", ">=", FC, "&", "p", ">", p, sep = " ")
    }
    con2 <- new_niden[(abs(new_niden$Fold_change) >= log2(FC) & new_niden$p.value >= -log10(p)),]
    if(dim(con2)[1] > 0){
      new_niden[(abs(new_niden$Fold_change) >= log2(FC) & new_niden$p.value >= -log10(p)),]$my_cor =
        paste("|FC|", ">=", FC, "&", "p", "<=", p, sep = " ")
    }
    if(dim(con2)[1] > 0){
      new_niden[(abs(new_niden$Fold_change) < log2(FC) & new_niden$p.value >= -log10(p)),]$my_cor =
        paste("|FC|", "<", FC, "&", "p", "<=", p, sep = " ")
    }

    ## plot
    myplot <- ggplot(data = new_niden, aes(x = Fold_change, y = p.value, label = Compoud_QueryID, col = my_cor)) +
      geom_point(alpha = 1, size = 1) +
      geom_text(hjust=0, vjust=0, alpha = 1) +
      scale_color_manual(values=c("steelblue1", "Tomato", "springgreen3", "Grey")) +
      theme_bw() +
      geom_vline(xintercept = log2(FC), col = "black", linetype = "dotted", size = 0.6) +
      geom_vline(xintercept = -log2(FC), col = "black", linetype = "dotted", size = 0.7) +
      geom_hline(yintercept = -log10(p), col = "black", linetype = "dotted", size = 0.6) +
      ggtitle("Not identified compound") +
      xlab("Log(2) fold change") +
      ylab("-Log(10) p-value") +
      labs(color = "")
  }


  #(4) final plots
  if(interactive == TRUE) {
    return(ggplotly(myplot))
  } else {
    return(myplot)
  }
}
