#' @title Element isotopes
#' @description check element isotope information
#' @author Yonghui Dong
#' @param S element, can be element symbol (i.e. C) or full name (i.e. Carbon).
#' Both Element symbol and full name are case insensitive.
#' @export
#' @examples
#'  E_iso('Na') # element symbol
#'  E_iso('nA') # element symbol, case insensitive
#'  E_iso('Carbon') # element full name
#'  E_iso('carBon') # element full name, case insensitive

# check isotopes
E_iso <- function(S) {
  #(1) read element data
  element <- as.data.frame(sysdata$element)
  element$Symbol <- as.character(element$Symbol)
  element$Name <- toupper(element$Name)
  element$Class <- toupper(element$Class)

  #(2) check isotopes
  S <- toupper(S)
  if (is.element(S, element$Class) == TRUE) {
    Symbol <- element$Symbol[element$Class == S]
    Abund <- element$Abund.[element$Class == S]
    Mass <- element$Mass[element$Class == S]
    isotable <- data.frame(element = Symbol, abundance = Abund, Mass = Mass)
    print(isotable)
  } else {
    if (is.element(toupper(S), element$Name) == TRUE){
      Symbol <- element$Symbol[element$Name == toupper(S)]
      Abund <- element$Abund.[element$Name == toupper(S)]
      Mass <- element$Mass[element$Name == toupper(S)]
      isotable <- data.frame(element = Symbol, abundance = Abund, Mass = Mass)
      print(isotable)
    } else {
      print('Not found. Use full name (i.e. Carbon) or symbol (i.e. C)')
    }
  }
}
