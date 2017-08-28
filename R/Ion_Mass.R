#' @title accurate ion mass
#' @description calculate accurate ion mass
#' @param m chemical formula of an ion
#' @param z charge, default value is 1
#' @param mode ionization mode, either positive '+' or negative '-'
#' @export
#' @examples
#'  mz('C7H7O4', z = 1, mode = '+')

# calculate accurate ion mass
mz <- function(m, z=1, mode = c('+', '-')) {
  ## read element data, and find the element with the highest abundance
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  e <- 0.000548597 # mass of an electron
  ## split the mass formula
  v1 <- strsplit(m, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  num <- as.numeric(v1[c(FALSE, TRUE)])
  if (length(atom) == length(num)) {
    options(digits=12)
    atom_mass <- element.max$Mass[match(atom, element.max$Class)]
    if (mode == '+') {
      accurate_mz <- (sum(atom_mass*num)-z*e)/z
      return(accurate_mz)
    } else if (mode == '-') {
      accurate_mz <- (sum(atom_mass*num)+z*e)/z
      return(accurate_mz)
    } else
      message ('Ion mode is not correct')
  }
  else
    message('Wrong chemical formula. If the numbers of some elements are missing?')
}
