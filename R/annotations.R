#' Find closest position in a reference
#'
#' For a list of positions, find the closest match in a given reference.
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
#' @param kallisto A \code{logical}. If set to TRUE, return a table that
#'  can be used as target_mapping for \code{sleuth_prep} or tx2gene for
#'  \code{tximport}.
#'
#' @importFrom biomaRt useMart getBM
#' @export
enst2symbol <- function(enst = NULL, version = F, kallisto = F){
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
  if(kallisto){
    e2g <- e2g[, c(ncol(e2g), 1), with = F]
    setnames(e2g, c("target_id", "GENE"))
  }
  return(e2g)
}


#' @importFrom biomaRt useMart getBM
ens_genes <- function(table, release = "hg38"){
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


#' GTEx color scheme
#'
#' This function maps tissues with their associated colour in the GTEx project
#'
#' @details
#' The colours are extracted from supplementary figure 5 of
#' "Genetic effects on gene expression across human tissues",
#' Nature volume 550, pages 204–213 (12 October 2017).
#'
#' @references  Consortium G. (2017). Genetic effects on gene expression across human tissues. Nature 550 204–213. 10.1038/nature24277
#'
#' @return A \code{data.table} mapping tissue, abbreviated name and colour.
#'
#' @export
GTEx_colors <- function(){
  colors <- c("#FF6600", "#FFAA00", "#33DD33", "#FF5555", "#FFAA99", "#FF0000",
              "#EEEE00", "#EEEE00", "#EEEE00", "#EEEE00", "#EEEE00", "#EEEE00",
              "#EEEE00", "#EEEE00", "#EEEE00", "#EEEE00", "#33CCCC", "#CC66FF",
              "#AAEEFF", "#EEBB77", "#CC9955", "#8B7355", "#552200", "#BB9988",
              "#9900FF", "#660099", "#AABB66", "#99FF00", "#AAAAFF", "#FFD700",
              "#FFAAFF", "#995522", "#AAFF99", "#DDDDDD", "#0000FF", "#7777FF",
              "#555522", "#778855", "#FFDD99", "#AAAAAA", "#006600", "#FF66FF",
              "#FF5599", "#FF00BB")
  abbreviations <- c("ADPSBQ", "ADPVSC", "ADRNLG", "ARTAORT", "ARTCRN", "ARTTBL",
                     "BRNACC", "BRNCDT", "BRNCHB", "BRNCHA", "BRNCTXA", "BRNCTXB",
                     "BRNHPP", "BRNHPT", "BRNNCC", "BRNPTM", "BREAST", "LCL", "FIBRBLS",
                     "CLNSGM", "CLNTRN", "ESPGEJ", "ESPMCS", "ESPMSL", "HRTAA", "HRTLV",
                     "LIVER", "LUNG", "MSCLSK", "NERVET", "OVARY", "PNCREAS", "PTTARY",
                     "PRSTTE", "SKINNS", "SKINS", "SNTTRM", "SPLEEN", "STMACH", "TESTIS",
                     "THYROID", "UTERUS", "VAGINA", "WHLBLD")
  tissues <- c("adipose - subcutaneous", "adipose - visceral (omentum)", "adrenal gland",
               "artery - aorta", "artery - coronary", "artery - tibial", "brain - anterior cingulate cortex (ba24)",
               "brain - caudate (basal ganglia)", "brain - cerebellar hemisphere [frozen]",
               "brain - cerebellum [paxgene]", "brain - cortex [paxgene]", "brain - frontal cortex (ba9) [frozen]",
               "brain - hippocampus", "brain - hypothalamus", "brain - nucleus accumbens (basal ganglia)",
               "brain - putamen (basal ganglia)", "breast - mammary tissue",
               "cells - ebv-transformed lymphocytes", "cells - transformed fibroblasts",
               "colon - sigmoid", "colon - transverse", "esophagus - gastroesophageal junction",
               "esophagus - mucosa", "esophagus - muscularis", "heart - atrial appendage",
               "heart - left ventricle", "liver", "lung", "muscle - skeletal",
               "nerve - tibial", "ovary", "pancreas", "pituitary", "prostate",
               "skin - not sun exposed (suprapubic)", "skin - sun exposed (lower leg)",
               "small intestine - terminal ileum", "spleen", "stomach", "testis",
               "thyroid", "uterus", "vagina", "whole blood")
  coltab <- data.table(tissue = tissues, abbreviation = abbreviations, color = colors)
  return(coltab)
}


