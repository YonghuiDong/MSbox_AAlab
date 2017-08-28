#' @title mass accuracy
#' @description calculate the mass accuracy of measured m/z
#' @param m measured m/z
#' @param t theoretical m/z
#' @export
#' @examples
#'  ppm(155.03383, 155.03388) # with m/z value
#'  ppm(155.03383, mz('C7H7O4', z = 1, mode = '+')) # with ion formula

# Calculate m/z accuracy
ppm <- function(m, t) {
    ppm <- (m-t)/t*10^6
    print(ppm)
}
