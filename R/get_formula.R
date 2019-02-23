#' @title Get compound formula and structure
#' @description get compound formula and structure from https://cactus.nci.nih.gov/chemical/structure
#' @author Yonghui Dong
#' @param chem, chemical name of the compound
#' @param representation, representation methods, formula is default
#' @param info, extra molecular information that users can query
#' @import xml2
#' @importFrom utils URLencode
#' @export
#' @examples
#' formula('malic acid')
#' formula(c('malic acid', 'citric acid', 'tartaric acid'))

formula <- function(chem, representation = 'formula', info = FALSE) {

  ##(1): display representation parameters
  if(info == TRUE) {message('More molecular information can be obtained by setting the representation parameter to:

[===========================================================================================]

(1) mw: molecular weight,
(2) monoisotopic_mass: monoisotopic_mass,
(3) h_bond_donor_count: number of hydrogen bond donors,
(4) h_bond_acceptor_count: number of hydrogen bond acceptors,
(5) h_bond_center_count: number of heavy atoms acting as hydrogen bond donor or acceptor,
(6) rule_of_5_violation_count: rule of 5 violation count,
(7) rotor_count,
(8) effective_rotor_count: Effective Rotor Count
(9) ring_count,
(10) ringsys_count,
(11) xlogp2,
(12) heteroatom_count,
(13) hydrogen_atom_count,
(14) heavy_atom_count,
(15) deprotonable_group_count,
(16) protonable_group_count,

[===========================================================================================]')}

  ##(2): query compound information
  root <- "https://cactus.nci.nih.gov/chemical/structure"
  ## define the lists
  url <- vector("list", length = length(chem))
  url_read <- vector("list", length = length(chem))
  url_result <- vector("list", length = length(chem))
  url_result2 <- vector("list", length = length(chem))

  for (i in 1:length(chem)) {
    ##(2.1) query compound formula
    url[[i]] <- paste(root, URLencode(chem[i]), representation, 'xml', sep = '/')
    url_read[[i]] <- read_xml(url[[i]])
    url_result[[i]] <- xml_text(xml_find_all(url_read[[i]], '//item'))
    ## check compound name
    if (identical(url_result[[i]], character(0)) == TRUE) {stop("Warning: compound name not found")}
    url_result2[[i]] <- url_result[[i]][[1]]
  }

  ## display compound other information
  names(url_result2) <- chem
  message('The ', representation, ' and compound structure(s) are as following:' )
  noquote(unlist(url_result2))
}
