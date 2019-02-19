#' fread file list
#'
#' fread a list of files and rbind them
#'
#' @param lf A \code{character} vector of filenames.
#' @param zcat A \code{logical}. Set to true if the files are zipped.
#' @param named A \code{logical}. Set to true if the filenames should be added
#'  as a column.
#' @param intersect A \code{logical}. Set to true if some elements have
#' different number of columns.
#'
#' @return A \code{data.table}.
#'
#' @export
freadl <- function(lf, zcat = F, named = T, intersect = T){
  fl <- vector('list', length(lf))
  for(i in seq_along(lf)){
    if(zcat)
      fri <- paste0("zcat ", lf[i])
    else
      fri <- lf[i]
    fl[[i]] <- fread(cmd = fri)
    if(named)
      fl[[i]][, filename := lf[i]]
  }
  if(intersect){
    nl <- Reduce(intersect, lapply(fl, names))
    ret <- rbindlist(lapply(fl, function(XX){XX[, nl, with = F]}))
  } else{
    ret <- rbindlist(fl)
  }
  return(ret)
}
