#' Find closest position in a reference
#'
#' @param pos A \code{numeric}. The list of position to match.
#' @param ref The \code{numeric} vector of postions to look into.
#' @export
find_nearest_pos <- function(pos, ref){
  idxs <- sapply(pos, function(XX) which.min(abs(XX - ref)))
  return(idxs)
}

#' Transcript to gene symbol
#'
#' Use biomaRt to obtain an annotation data.table that maps from ensembl
#' transcript to gene symbol.
#'
#' @param enst A \code{character} vector or NULL. The annotation will be
#' subset for the trancripts specified.
#' @param version A \code{logical}. If set to TRUE, return transcript version.
#'
#' @importFrom biomaRt useMart getBM
#' @export
enst2symbol <- function(enst = NULL, version = F){
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  if(version){
    atts <- c("hgnc_symbol", "ensembl_transcript_id", "ensembl_transcript_id_version")
  } else{
    atts <- c("hgnc_symbol", "ensembl_transcript_id")
  }
  e2g <- data.table(getBM(attributes = atts, mart = mart))
  setnames(e2g, c("hgnc_symbol", "ensembl_transcript_id"), c("GENE", "transcript"))
  if(!is.null(enst)){
    e2g <- e2g[transcript %in% enst]
  }
  return(e2g)
}
