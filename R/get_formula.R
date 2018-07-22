#' @title Get compound formula and structure
#' @description get compound formula and structure from https://cactus.nci.nih.gov/chemical/structure
#' @param chem, chemical name of the compound
#' @param representation, representation methods, formula is default
#' @param info, extra molecular information that users can query
#' @import xml2
#' @import magick
#' @importFrom utils URLencode
#' @export
#' @examples
#' formula('malic acid')

formula <- function(chem, representation = 'formula', info = FALSE) {
  if(info == TRUE) {message('more molecular information can be obtained by setting the
representation parameter to:mw, monoisotopic_mass, h_bond_donor_count, h_bond_acceptor_count,
h_bond_center_count, rule_of_5_violation_count, rotor_count, effective_rotor_count, ring_count,
ringsys_count, xlogp2, heteroatom_count, hydrogen_atom_count, heavy_atom_count,deprotonable_group_count, protonable_group_count')}

  root <- "https://cactus.nci.nih.gov/chemical/structure"
  url <- paste(root, URLencode(chem), representation, 'xml', sep = '/')
  ## restrict query time
  Sys.sleep(1.1)
  myread <- read_xml(url)
  myresult <- xml_text(xml_find_all(myread, '//item'))
  ## check compound name
  if (identical(myresult, character(0)) == TRUE) {stop("Warning: compound name not found")}
  myresult2 <- myresult[[1]]
  url_image <- paste(root, URLencode(chem), 'image', sep = '/')
  chem_image <- print(image_read(url_image), info = F)
  return(myresult2)
  return(chem_image)
}


