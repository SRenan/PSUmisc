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


#'
#' @details
#' This function is for the human genome only.
#'
#' @importFrom biomaRt useMart getBM
pos2gene <- function(table, release = "hg19"){
  release <- tolower(release)
  if(release %in% c("grch38", "hg38", "latest")){
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  } else if(release %in% c("grch37", "hg19")){
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  } else if(release %in% c("hg17")){

  } else{
    stop("Incorrect release number")
  }
  atts <- c("hgnc_symbol", "ensembl_transcript_id", "chromosome_name", "start_position", "end_position")
  p2g <- data.table(getBM(attributes = atts, mart = mart))
  setnames(p2g, c("hgnc_symbol", "start_position", "end_position"))
}

pos <- c(54795432, 131511098)

# Inputs: A table of positions + A table of genes with their ranges
# Output: A table with all SNPs that were found within a gene

#' Find genes overlapping a position
#'
#' @param target a \code{data.table} with two columns, ID and positions to match.
#' @param interval a \code{data.table} with at least three columns "GENE",
#'  "start" and "end".
#'
#' @note
#' When there are multiple matches in the \code{interval}, this function will
#' return all genes overlapping the query.
#'
#' @return A \code{data.table} with all positions that were within the given
#'  intervals along with the gene that covers them and its coordinates
#'
#'
#' @export
pos2gene <- function(target, interval){
  pos <- copy(target)
  gene <- copy(interval)
  # TODO:
  # a) Handle chromosome
  # b) Fetch genes from biomaRt

  # INPUTS
  # Assume first column is an ID and second is the positions
  setnames(pos, c("ID", "POS"))
  pos[, start := POS]
  pos[, end := POS]
  setkey(pos, start, end)
  # Assume genes have columns GENE, start, end
  olps <- foverlaps(gene, pos, type = "any")
  ret <- olps[!is.na(ID), list(ID, POS, GENE, i.start, i.end)]
  setnames(ret, c("i.start", "i.end"), c("start", "end"))
  return(ret)
}


