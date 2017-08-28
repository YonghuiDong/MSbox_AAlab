#' @title Element isotopes
#' @description Check element isotope information
#' @param S element, can be element symbol (i.e. C) or full name (i.e. Carbon).
#' Element symbol is case sensitive, which full name is case insensitive.
#' @export
#' @examples
#'  E_iso('C') # element symbol
#'  E_iso('Carbon') # element full name
#'  E_iso('carBon') # element full name

# check isotopes
E_iso <- function(S) {
  # read element data
  element <- as.data.frame(sysdata$element)
  element$Symbol <- as.character(element$Symbol)
  element$Name <- toupper(element$Name)
  # check isotopes
  if (is.element(S, element$Class) == TRUE) {
    Symbol <- element$Symbol[element$Class == S]
    Abund <- element$Abund.[element$Class == S]
    Mass <- element$Mass[element$Class == S]
    isotable <- rbind(element = Symbol, abundance = Abund, Mass = Mass)
    print(isotable)
  } else {
    if (is.element(toupper(S), element$Name) == TRUE){
      Symbol <- element$Symbol[element$Name == toupper(S)]
      Abund <- element$Abund.[element$Name == toupper(S)]
      Mass <- element$Mass[element$Name == toupper(S)]
      isotable <- rbind(element = Symbol, abundance = Abund, Mass = Mass)
      print(isotable)
    } else {
      print('Not found. Use full name (i.e. Carbon) or symbol (i.e. C)')
    }
  }
}
