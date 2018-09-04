#'
#' fread a list of files from
#'
#' @param lf A \code{character} vector of filenames.
#' @param zcat A \code{logical}. Set to true if the files are zipped.
#' @param named A \code{logical}. Set to true if the filenames should be added
#'  as a column.
#'
#' @export
freadl <- function(lf, zcat = F, named = T){
  fl <- vector('list', length(lf))
  for(i in seq_along(lf)){
    if(zcat)
      fri <- paste0("zcat ", lf[i])
    else
      fri <- lf[i]
    fl[[i]] <- fread(fri)
    if(named)
      fl[[i]][, filename := lf[i]]
  }
  ret <- rbindlist(fl)
}
