#' @title Get compound smiles
#' @description get compound formula and structure from https://cactus.nci.nih.gov/chemical/structure
#' @author Yonghui Dong
#' @param chem, chemical name of the compound
#' @import xml2
#' @importFrom utils URLencode
#' @export
#' @examples
#' get_smiles('malic acid')
#' get_smiles(c('malic acid', 'citric acid', 'tartaric acid'))

get_smiles <- function(chem) {

  ##(1): query compound information
  root <- "https://cactus.nci.nih.gov/chemical/structure"
  ## define the lists
  url <- vector("list", length = length(chem))
  url_read <- vector("list", length = length(chem))
  url_result <- vector("list", length = length(chem))
  url_result2 <- vector("list", length = length(chem))
  missing <- rep(NA, length(chem)) # count the unassigned number

  for (i in 1:length(chem)) {
    ##(2.1) query compound formula
    url[[i]] <- paste(root, URLencode(chem[i]), smiles, 'xml', sep = '/')
    url_read[[i]] <- read_xml(url[[i]])
    url_result[[i]] <- xml_text(xml_find_all(url_read[[i]], '//item'))
    ## check compound name
    if (identical(url_result[[i]], character(0)) == TRUE) {
      url_result2[[i]] = "unknown"
      missing[i] = i # record the missing index
      } else {
        url_result2[[i]] <- url_result[[i]][[1]]
      }
  }

  ## display compound other information
  names(url_result2) <- chem
  missing2 <- missing[!is.na(missing)]
  message('The ', representation, ' are as following:' )
  if (length(missing2) > 0) {
    message("Attention: smiles of ", length(missing2), " compound(s) ", "are not assigned")
  }
  noquote(unlist(url_result2))
}





