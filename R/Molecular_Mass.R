#' @title molecular mass
#' @description calculate accurate molecular mass
#' @param F chemical formula, case insensitive
#' @importFrom stats aggregate
#' @export
#' @examples
#'  M_mass('C7H6O4')
#'  M_mass('c7H6O4') # case insensitive

M_mass <- function(F) {
  #(1) read element data, and find the element with the highest abundance
  options(digits = 12)
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  element.max$Class <- toupper(element.max$Class)
  #(2) split the mass formula
  v1 <- strsplit(F, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  ## convert the first letter of atom to upper case. case insensitive.
  atom <- toupper(atom)
  ## check the input
  atom_match <- match(atom, element.max$Class)
  if(any(is.na(atom_match))) {stop("Warning: certain element not found")}
  num <- as.numeric(v1[c(FALSE, TRUE)])
  if (length(atom) == length(num)) {
    atom_mass <- element.max$Mass[match(atom, element.max$Class)]
    accurate_mass <- sum(atom_mass*num)
    ## set 6 decimal
    accurate_mass = formatC(accurate_mass, digits = 6, format = "f")
    accurate_mass <- as.numeric(accurate_mass)
    return(accurate_mass)
  }
  else
    message('Wrong chemical formula. Numbers of some elements missing?')
}
