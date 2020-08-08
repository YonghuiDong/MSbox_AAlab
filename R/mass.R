#' @title molecular mass
#' @description calculate accurate molecular mass
#' @author Yonghui Dong
#' @param F chemical formula, case insensitive
#' @importFrom stats aggregate
#' @importFrom crayon cyan bgMagenta
#' @export
#' @examples
#'  mass('C7H6O4')
#'  mass('c7H6O4') # case insensitive
#'  mass(c('K1', 'C5H8', 'nA20')) # vector input

mass <- function(F) {

  #(1) read element data, and find the element with the highest abundance
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  element.max$Class <- toupper(element.max$Class)

  #(2) allow the function to a vector of input
  if (length(F) > 1) {
    F = as.list(F)
    acc_mass <- sapply(F, mass)
    names(acc_mass) <- F
    return(acc_mass)
  }

  #(3) main function
  ## split the mass formula
  v1 <- strsplit(F, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  ## convert the first letter of atom to upper case. case insensitive.
  atom <- toupper(atom)
  ## check the input
  ###(1) check element
  atom_match <- match(atom, element.max$Class)
  if(any(is.na(atom_match))) {stop("Warning: certain element not found")}
  ## special elements that we need to pay attention
  ## for example, we want to know mass of C(carbon)O(oxygen), if your forget to write as C1O1,
  ## it will be calculated for element Co (cobalt)
  Exlist <- c("Co", "Cs", "Cu", "Ho", "In", "Ni", "Os", "Sn", "Si")
  exmatch <- Exlist[which(toupper(Exlist) %in% atom)]
  if(length(exmatch) > 0)
  {cat(cyan("Attention: are following elements indeed in your formula: "), bgMagenta(exmatch), "\n")}
  ###(2) check element number
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
