#' @title accurate ion mass
#' @description calculate accurate ion mass
#' @param m chemical formula of an ion, case insensitive
#' @param z charge
#' @export
#' @examples
#'  mz('C7H7O4', z = 1)
#'  mz('C10H6Cl1', z = -1)
#'  mz('C7h7O4', z = 1) # case insensitive

# calculate accurate ion mass
mz <- function(m, z) {
  # (1) issue warnings
  if(z == 0) {stop("Warning: charge z = 0 ?")}
  if(z%%1 != 0) {stop("Warning: charge z must be integer")}
  if(is.numeric(z) == FALSE) {stop("Warning: charge z shoule be numeric!")}
  ## read element data, and find the element with the highest abundance
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  e <- 0.000548597 # mass of an electron
  ## split the mass formula
  v1 <- strsplit(m, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  # convert the first letter of atom to upper case. case insensitive.
  atom <- paste(toupper(substr(atom, 1, 1)), substr(atom, 2, nchar(atom)), sep="")
  num <- as.numeric(v1[c(FALSE, TRUE)])
  if (length(atom) == length(num)) {
    options(digits=12)
    atom_mass <- element.max$Mass[match(atom, element.max$Class)]
    accurate_mz <- (sum(atom_mass*num)-z*e)/abs(z)
    return(accurate_mz)
  }
  else
    message('Wrong chemical formula. Are numbers of some elements missing?')
}
