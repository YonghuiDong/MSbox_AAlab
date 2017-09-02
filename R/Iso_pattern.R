#' @title Isotope pattern distribution
#' @description Calculate isotope pattern distribution of given mass formula(e).
#' @note This script is adapted from R package 'enviPat', all the credits go to the original author
#' @param F chemical formula, must be written in capital letters.
#' @param threshold only display the m/z values with relative abundance over the threshold,
#'  default value = 0.01.
#' @param z charge.
#' @param rel_to 0: relative to highest m/z; 1: relative to monoisotopic peak.
#' @param return_iso_calc_amount ???
#' @export
#' @examples
#'  iso_pattern('C7H7O4',z = 1)
#'  iso_pattern(c('C7H7O4', 'C4H7O6'), z = c(1, -1))


#'  library(plotly)
#'  out = iso_pattern('C700H7O4K30',z = 1)
#'  out = as.data.frame(out$C700H7O4K30)
#'  plot_ly(out) %>% add_trace(x = out$`m/z`, y = out$abundance, type = 'bar', width = .001)

iso_pattern <-
  function(
    F,
    threshold = 0.01,
    z,
    rel_to = 0,
    return_iso_calc_amount=FALSE
  ){

    ## (1) issue warnings
    if(threshold>100 | threshold<0){stop("WARNING: invalid threshold; 0<=threshold<100.\n")}
    if((length(z)!=length(F)) & length(z)>1){stop("length of charge does not match number of chemforms!\n")}
    if(any(z==0)) {stop("WARNING: charge=0?")}
    if(any(is.numeric(z)==FALSE)) {stop("WARNING: charge should be numeric!")}
    if(length(z)==1 & length(F)>1){z <- rep(z,length(F))}
    options(digits=10);
    if(!any(rel_to==c(0,1,2,3,4))){stop("invalid rel_to")}
    if(return_iso_calc_amount=="TRUE"){return_iso_calc_amount2=1}else{return_iso_calc_amount2=0}

    ## (2) run isotope pattern generator
    e <- 0.000548597 # mass of an electron
    load('R/sysdata.rda')
    isotopes <- as.data.frame(sysdata$isopattern)
    pattern<-list(0)
    for(i in 1:length(F)){
        out <- .Call( "iso_pattern",
                      s1 = as.character(F[i]),   # chemical formula
                      pl = as.integer(1E6),             # number of peaks to be reserved for
                      t1 = as.double(threshold),        # relative intensity cutoff
                      iso_list_elem = as.character(isotopes[,1]),  # isotope list: Element
                      iso_list_iso = as.character(isotopes[,2]),  # isotope list: Isotope
                      iso_list_mass = as.numeric(isotopes[,3]),  # isotope list: Isotope mass
                      iso_list_abu = as.numeric(isotopes[,4]),  # isotope list: Isotope abundance
                      rtm = as.integer(rel_to),     # 0:relative to highest, 1:relative to mono peak
                      rica = as.integer(return_iso_calc_amount2)
        )
      # parse output ###########################################################
      if(length(out[[1]])==0){
        pattern[[i]]<-"error";
      }else{
        if(return_iso_calc_amount2){
          pattern[[i]]<-out
        }else{
          out2<-out[order(out[,1],decreasing=FALSE),,drop=FALSE]
          colnames(out2)[1]<-"m/z"
          out2[,1]<-(out2[,1]-(z[i]*e))/abs(z[i])  # electrone mass
          pattern[[i]]<-out2

          plot(out2[,1],out2[,2],type="h",
               xlab="m/z",ylab="Relative abundance",main=names(pattern)[i])

        }
      }
    }
    names(pattern)<-as.character(F);

    ## (3) output
    return(pattern)
  }
