#' @title Common adducts
#' @description calculate common adduct ions in pos or neg ion mode
#' @param F chemical formula, must be written in capital letters
#' @param mode ionization mode, either positive '+' or negative '-'
#' @export
#' @examples
#' adduct('C1H4',mode = '-')
#' adduct('C1H4',mode = '+')

adduct <- function(F, mode = c('+', '-')){
  if (mode == '+') {
    pos_ion <- c('Li1', 'H3O1', 'N1H4', 'Na1', 'K1', 'C1H5O1', 'C2H4N1')

    M_adduct <- vector(mode="character", length=length(pos_ion))
    for (i in 1: length(pos_ion)) {
      M_adduct[i] = mz(paste(F, pos_ion[i], sep = ''), z = 1, mode = '+')
    }
    adduct_info <- cbind(adduct = c('M+Li','M+H3O', 'M+NH4', 'M+Na', 'M+K', 'M+H+CH3OH',
                                    'M+H+CH3CN'),
                         mz = M_adduct,
                         source = c('Lithium salts', 'Ammonia/NH4OH', 'Water/Acids',
                                    'Sodium salts', 'Potassium salts', 'Methanol', 'Acetonitrile')
    )
    return(adduct_info)
  } else if (mode == '-') {
    neg_ion <- c('C1H3O1', 'C1O2H1', 'C2H3O2','C2F3O2', 'Cl1')
    M_adduct <- vector(mode="character", length=length(neg_ion))
    for (i in 1: length(neg_ion)) {
      M_adduct[i] = mz(paste(F, neg_ion[i], sep = ''), z = 1, mode = '-')
    }
    M_adduct <- c(M_adduct, M_mass(F) + 36.965903 + 0.000548597)
    adduct_info <- cbind(adduct = c('M-H+CH3OH','M-H+HCO2H','M-H+CH3CO2H', 'M-H+CF3CO2H',
                                    'M+Cl(35)', 'M+Cl(37)'),
                         mz = M_adduct,
                         source = c('Methanol', 'Formic acid', 'Acetic acid', 'TFA',
                                    'Chlorinated solvent', 'Chlorinated solvent')
                         )
    return(adduct_info)
  } else
    message ('Ion mode is not correct')
}
