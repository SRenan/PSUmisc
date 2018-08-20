#'
#' fread a list of files from
#'
#' @export
freadl <- function(lf, zcat = F){
  fl <- vector('list', length(lf))
  for(i in seq_along(lf)){
    if(zcat)
      fri <- paste0("zcat ", lf[i])
    else
      fri <- lf[i]
    fl[[i]] <- fread(fri)
  }
  ret <- rbindlist(fl)
}
