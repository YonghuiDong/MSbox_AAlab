#' @title fold change
#' @description calculate fold change among different samples.
#' @param x dataframe
#' @importFrom utils combn
#' @return a dataframe with mean values and fold changes
#' @export
#' @examples
#' vars <- 1000
#' samples <- 50
#' Groups <- 3
#' dat <- replicate(vars, runif(n = samples))
#' f <- rep_len(1:Groups, samples)
#' f <- LETTERS[f]
#' dat <- data.frame(dat, Group = f)
#' ret <- fold(dat)

 fold <- function(x){
  f <- x$Group
  i <- split(1:nrow(x), f)
  ## rm Group info in order to calculate colMeans
  x$Group = NULL
  mean_int <- sapply(i, function(i){colMeans(x[i, ])})
  x <- t(mean_int)
  j <- combn(levels(f), 2)
  f_change1 <- x[j[1,],] / x[j[2,],]
  f_change2 <- x[j[2,],] / x[j[1,],]
  ## remove NaN in f_change Matrix
  f_change <- rbind(f_change1, f_change2)
  f_change[is.nan(f_change)] <- 0
  rownames(f_change) <- c(paste("Fold_", j[1,], "_VS_", j[2,], sep = ''), paste("Fold_", j[2,], "_VS_", j[1,], sep = '') )
  ret <- as.data.frame(t(f_change))
  return(ret)
}
