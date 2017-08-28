#' @title molecular mass
#' @description calculate accurate molecular mass
#' @param F chemical formula, must be written in capital letters
#' @export
#' @examples
#'  M_mass('C7H6O4')

# calculate accurate molecular mass
M_mass <- function(F) {
  ## read element data, and find the element with the highest abundance
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  ## split the mass formula
  v1 <- strsplit(F, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  num <- as.numeric(v1[c(FALSE, TRUE)])
  if (length(atom) == length(num)) {
    options(digits=12)
    atom_mass <- element.max$Mass[match(atom, element.max$Class)]
    accurate_mass <- sum(atom_mass*num)
    return(accurate_mass)
  }
  else
    message('Wrong chemical formula. Numbers of some elements missing?')
}
