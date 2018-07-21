#' @title Get compound formula
#' @description get compound formula from https://cactus.nci.nih.gov/chemical/structure
#' @param chem, chemical name of the compound
#' @param representation, representation methods, formula is default
#' @import xml2
#' @importFrom utils URLencode
#' @export
#' @examples
#' formula('malic acid')

formula <- function(chem, representation = 'formula') {
  ref_url <- "https://cactus.nci.nih.gov/chemical/structure"
  chem_url <- paste(ref_url, URLencode(chem), representation, 'xml', sep = '/')
  ## restrict query time
  Sys.sleep(1.1)
  myread <- read_xml(chem_url)
  ## check compound name
  if (identical(myread, character(0)) == TRUE) {stop("Warning: compound name not found")}
  myresult <- xml_text(xml_find_all(myread, '//item'))[[1]]
  return(myresult)
}
