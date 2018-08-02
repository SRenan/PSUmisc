#' Process kallisto outputs
#'
#' @param bstrapdir A \code{character}. The path of a kallisto bootstrap directory.
#'
#' @details
#' This only reads the tsv files to get raw counts.
#'
#' @export
readKali <- function(bstrapdir){
  samples <- list.files(bstrapdir)
  abfiles <- file.path(bstrapdir, samples, "abundance.tsv")
  exl <- file.exists(abfiles)
  if(!all(exl)){
    msg <- paste("There are", sum(!exl), "missing abundance files")
    warning(msg)
  }
  samples <- samples[exl]
  paths <- file.path(bstrapdir, samples)
  abunds <- file.path(paths, "abundance.tsv")
  ablist <- vector("list", length(abunds))
  for(i in seq_along(abunds)){
    ablist[[i]] <- fread(abunds[i])
    ablist[[i]][, sample := samples[[i]]]
  }
  ablist <- rbindlist(ablist)
  tpms <- merge(ablist, map, by = "sample")
  tpms[, transcript := gsub("\\..*", "", target_id)]

  return(tpms)
}
