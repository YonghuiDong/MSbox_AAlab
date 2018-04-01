#' @title mass accuracy
#' @description calculate the mass accuracy of measured m/z. lazy input allowed.
#' @param m measured m/z
#' @param t theoretical m/z
#' @param lazy if lazy input is allowed
#' @export
#' @examples
#'  ppm(155.03383, 155.03388) # with m/z value
#'  ppm(155.03383, .03388) # lazy input when the integer parts of m and t are the same
#'  ppm(155.03384, mz('C7H7O4', z = 1)) # with ion formula

ppm <- function(m, t, lazy = TRUE) {
  # To disable scientific notation. This is important, otherwise it will bring bug when during spliting
  options(scipen=999)
  # split the theorectical m/z t, and extract the integer part
  t1 = unlist(strsplit(as.character(t),"\\."))[1]
  t1 <- as.numeric(t1)
  # split the measured m/z m, and extract the integer part
  m1 = unlist(strsplit(as.character(m),"\\."))[1]
  m1 <- as.numeric(m1)
  # split measured m/z m, and extract the fractional part
  m2 = unlist(strsplit(as.character(m),"\\."))[2]
  if(is.na(m2) == TRUE){m2 = 0}
  m2 <- paste0(".", m2)
  m2 <- as.numeric(m2)
  # split the theorectical m/z t, and extract the fractional part
  t2 = unlist(strsplit(as.character(t),"\\."))[2]
  if(is.na(t2) == TRUE){t2 = 0}
  t2 <- paste0(".", t2)
  t2 <- as.numeric(t2)
  # if the integer part of t is 0, which means it has the same inerger part as m
  if (t1 == 0 & lazy == TRUE) {
    mz_dif <- m2 - t2
    mz_dif <- formatC(mz_dif, digits = 5, format = "f")
    ppm <- (m2-t2)/(t2+m1)*10^6
    ppm <- formatC(ppm, digits = 4, format = "f")
    result <- as.data.frame(cbind(mz_dif = mz_dif, ppm = ppm))
    print.data.frame(result)
  } else {
    mz_dif <- m - t
    mz_dif <- formatC(mz_dif, digits = 5, format = "f")
    ppm <- (m-t)/t*10^6
    ppm <- formatC(ppm, digits = 4, format = "f")
    result <- as.data.frame(cbind(mz_dif = mz_dif, ppm = ppm))
    print.data.frame(result)
  }
}
