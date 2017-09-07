#' @title Isotope labelled molecular mass
#' @description Calculate isotope labelled molecular mass
#' @param F, chemical formula, case insensitive
#' @param iso, labelled elements, case insensitive
#' @param z charge
#' @importFrom stats aggregate
#' @export
#' @examples
#' Iso_mz(F = 'C7H6O4', iso = '[13]C2[2]H3', z = -1) # Two 13C and three 2H are labled

Iso_mz <- function(F, iso, z) {
  element <- as.data.frame(sysdata$element)
  # replace '()' with '[]'
  element$Symbol <- chartr("()", "[]", element$Symbol)
  # change the formact i.e. C[13] to [13]C
  element$Symbol <- gsub("^(.+)(\\[[0-9]+\\])$", "\\2\\1", element$Symbol)
  element$Abund.<- as.numeric(element$Abund.)
  element.agg <- aggregate(Abund. ~ Class, element, max)
  element.max <- merge(element.agg, element)
  # split iso
  grx <- gregexpr("\\[.+?\\].+[[:digit:]]?",  iso)
  let <- do.call(c, regmatches(iso, grx))
  grx <- gregexpr("\\[.+?\\].+([[:digit:]]+)",  iso)
  out <- do.call(c, regmatches(iso, grx))
  num <- gsub(".+\\][[:alpha:]]+", "", out)
  num <- as.numeric(num)
  # match iso
  iso_atom_mass <- element$Mass[match(let, element$Symbol)]
  # match the monoisotopic atom
  iso_atom_class <- element$Class[match(let, element$Symbol)]
  atom_mass <- element.max$Mass[match(iso_atom_class, element.max$Class)]
  # calculate mass difference
  mass_dif <- (iso_atom_mass - atom_mass)
  # calculate total increased mass
  mass_total_increase <- sum(mass_dif * num)

  # calculate isto_mass
  iso_mz <- mz(F, z) + mass_total_increase

  return(iso_mz)
}
