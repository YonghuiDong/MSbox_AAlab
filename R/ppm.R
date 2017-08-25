#' @title mass accuracy
#' @description calculate the mass accuracy of measured m/z
#' @author Yonghui Dong
#' @param m measured m/z
#' @param t theoretical m/z
#' @examples
#' ppm(155.03383, 155.03388)

# Calculate m/z accuracy
ppm <- function(m, t) {
    ppm <- (m-t)/t*10^6
    print(ppm)
}
