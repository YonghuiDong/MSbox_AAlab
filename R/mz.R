#' @title accurate ion mass
#' @description calculate accurate ion mass
#' @author Yonghui Dong
#' @param m chemical formula of an ion, case insensitive
#' @param z charge
#' @importFrom crayon cyan bgMagenta
#' @export
#' @examples
#'  mz('C7H7O4', z = 1)
#'  mz('C10H6Cl1', z = -1)
#'  mz('C7h7O4', z = 1) # case insensitive
#'  mz(c('C7H7O4', 'c1'), z = -1) # vector input

mz <- function(m, z) {
  options(digits = 12)

  #(1) check input
  if(z == 0) {stop("Warning: charge z = 0 ?")}
  if(z%%1 != 0) {stop("Warning: charge z must be integer")}
  if(is.numeric(z) == FALSE) {stop("Warning: charge z shoule be numeric!")}

  #(2) read element data, and find the element with the highest abundance
  element <- as.data.frame(sysdata$element)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  element.max$Class <- toupper(element.max$Class)
  e <- 0.000548597 # mass of an electron

  #(3) allow the function to a vector of input
  if (length(m) > 1) {
    m = as.list(m)
    acc_mz <- sapply(m, mz, z = z)
    names(acc_mz) <- m
    return(acc_mz)
  }

  #(4) main function
  ## format input
  ## split the mass formula
  v1 <- strsplit(m, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  # convert the first letter of atom to upper case. case insensitive.
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
    accurate_mz <- (sum(atom_mass*num)-z*e)/abs(z)
    accurate_mz = formatC(accurate_mz, digits = 6, format = "f")
    accurate_mz <- as.numeric(accurate_mz)
    return(accurate_mz)
  }
  else
    message('Wrong chemical formula. Are numbers of some elements missing?')
}
