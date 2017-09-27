#' @title Isotope labelled molecular mass
#' @description Calculate isotope labelled molecular mass
#' @param F, chemical formula, case insensitive
#' @param iso, labelled elements, case insensitive
#' @importFrom stats aggregate
#' @importFrom stringr str_match_all
#' @export
#' @examples
#' Iso_mass(F = 'C7H6O4', iso = '[13]C2[2]H3') # Two 13C and three 2H are labled

Iso_mass <- function(F, iso) {
  element <- as.data.frame(sysdata$element)
  #(1) format the database
  ## replace '()' with '[]'
  element$Symbol <- chartr("()", "[]", element$Symbol)
  ## change the formact i.e. C[13] to [13]C
  element$Symbol <- gsub("^(.+)(\\[[0-9]+\\])$", "\\2\\1", element$Symbol)
  element$Symbol <- toupper(element$Symbol)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  #(2) split iso
  grx <- str_match_all(iso, "(\\[\\d+\\]\\p{L}+)(\\d+)")
  ## into letter
  let <- grx[[1]][,2]
  let <- toupper(let)
  ## into number
  iso_num <- grx[[1]][,3]
  iso_num <- as.numeric(iso_num)
  #(3) check the number of labelled element should not exceed the number of that element
  ## prepare iso infomration
  iso_atom_class <- as.character(element$Class[match(let, element$Symbol)])
  iso_infor <- data.frame(iso_atom_class, iso_num)
  ## prepare molecule information
  ## split the mass formula
  v1 <- strsplit(F, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)[[1]]
  atom <- v1[c(TRUE, FALSE)]
  ## convert the first letter of atom to upper case. case insensitive.
  atom <- paste(toupper(substr(atom, 1, 1)), substr(atom, 2, nchar(atom)), sep="")
  total_num <- as.numeric(v1[c(FALSE, TRUE)])
  F_infor <- data.frame(atom, total_num)
  ## match iso in molecule
  index <- match(iso_infor$iso_atom_class, F_infor$atom)
  if(any(is.na(index)) == TRUE) {stop ("Warning: certain labelled element not found in the molecule")}
  if(any(iso_infor$iso_num > F_infor$total_num[index]) == TRUE)
    {stop ("Warning: The number of certain labelled element exceeds the max number in the molecule")}
  #(4) Calculate the iso mass
  # match iso
  iso_atom_mass <- element$Mass[match(let, element$Symbol)]
  # match the monoisotopic atom
  iso_atom_class <- element$Class[match(let, element$Symbol)]
  atom_mass <- element.max$Mass[match(iso_atom_class, element.max$Class)]
  # calculate mass difference
  mass_dif <- (iso_atom_mass - atom_mass)
  # calculate total increased mass
  mass_total_increase <- sum(mass_dif * iso_num)

  # calculate isto_mass
  iso_mass <- M_mass(F) + mass_total_increase

  return(iso_mass)
}
