library(gtools)
library(ggplot2)
library(data.table)

#' Draw manhattan plot
#'
#' @param chr Vector of chromosome names
#' @param pos Vector of genomic positions
#' @param pval Vector of P-values from association test
#'
#' @importFrom ggplot2 ggplot aes theme theme_bw element_blank
#' @importFrom ggplot2 geom_point geom_hline geom_rect
#' @importFrom ggplot2 scale_fill_manual scale_colour_manual scale_x_continuous
#' @importFrom ggplot2 xlab ylab
#' @importFrom gtools mixedorder
#' @export
ggman <- function(chr, pos, pval, names = NULL, minpval = 1e-50, signif = 5e-8, lowp = 5e-3){
  chr <- as.character(chr)
  pos <- as.numeric(pos)

  nchr <- length(unique(chr))
  mantab <- data.table(chr, pos, pval)
  mantab[, posmin := min(pos), by = "chr"]
  mantab[, posmax := max(pos), by = "chr"]

  chrtab <- unique(mantab[, list(chr, posmin, posmax)])[mixedorder(chr)]
  chrtab[, posshift := c(0, cumsum(posmax))[1:nchr]]
  chrtab[, poslabel := (posmax+posmin)/2 + posshift]
  chrtab[, chrcolor := rep_len(c("cyan", "gray80"), length.out = nchr)]

  tab <- merge(mantab[, list(chr, pos, pval)], chrtab, by = "chr")
  tab[pval < minpval, pval := minpval]
  tab[, wgpos := pos + posshift]
  tab[, log10p := -log10(pval)]

  scm <- scale_colour_manual(values = chrtab$chrcolor, guide = "none")
  sxc <- scale_x_continuous(breaks = chrtab$poslabel, labels = chrtab$chr)
  theme_man <- theme(panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank())
  p <- ggplot(tab[pval < lowp]) + geom_point(aes(x=wgpos, y=log10p, color = chrcolor)) +
    scm + sxc + theme_bw() + theme_man +
    xlab("Chromosome") + ylab(expression(-log[10](p-value)))
  if(signif > 0)
    p <- p + geom_hline(yintercept = -log10(signif), linetype = 2)
  if(lowp < 1){
    sfm <- scale_fill_manual(values = chrtab$chrcolor, guide = "none")
    p <- p + geom_rect(data = chrtab, aes(xmin = posshift, xmax = posshift + posmax, ymin = 0, ymax = -log10(lowp), fill = chrcolor)) +
      sfm
  }
  print(p)
  return(p)
}

#' Draw QQ-plot
#'
#' @importFrom ggplot2 geom_abline
#' @export
ggqqp <- function(pval, alpha = 1e-2){
  # The pvalues are assumed to follow a chisq with 1df
  lout <- length(pval)
  obs <- sort(pval)
  sigobs <- obs[obs < alpha]
  nsig <- length(sigobs)
  nunsig <- lout - nsig

  idx <- c(1:5, nunsig:lout)

  # exp <- seq(0, 1, length.out = lout)
  exp <- seq(1, lout)
  exp <- (exp - 0.5) / lout


  lambda <- round(qchisq(median(obs), df= 1, lower.tail = F)/qchisq(0.5, df=1, lower.tail = F), 3)
  lambda_dt <- data.table(lambda)
  dt <- data.table(obs, exp)
  dt <- dt[order(exp, decreasing = T)]
  p <- ggplot(dt[idx], aes(x = -log10(exp), y = -log10(obs))) + geom_point() + geom_abline() +
         geom_label(data = lambda_dt, aes(label = paste("lambda ==", lambda), x = 0, y = 0), parse = T, hjust = 0, vjust = 0)
  p <- p + theme_bw()
  print(p)
  return(p)
}


#' Isolate most significant SNP at each locus
#' @export
topSNP <- function(pos, pval, ret = "idx", locus_len = 1e6){
  # Order positions
  if(any(sort(pos) != pos)){
    stop("The positions are not ordered")
  }
  if(length(pos) == 1){
    return(1)
  }
  # Split the data into loci
  loci_start <- c(1, which(diff(pos) > locus_len)+1)
  loci_end <- c(loci_start - 1, length(pos))[-1]
  # Find max loci between start pos and end
  idx <- sapply(seq_along(loci_start), function(XX){which.min(pval[loci_start[XX]:loci_end[XX]])})
  idx <- idx + loci_start - 1
  if(ret == "idx"){
    ret <- idx
  } else if(ret == "pos"){
    ret <- pos[idx]
  }
  return(ret)
}

#' Subset a table for the top SNP in each locus
#'
#' @param table A \code{data.table} of significantly associated SNPs, with the
#'  columns listed in the next 3 arguments. Additional columns are preserved
#' @param pos A \code{character}. The name of the column of positions in \code{table}.
#' @param pval A \code{character}. The name of the column of p-values in \code{table}.
#' @param group A \code{character}. The name of the column used for grouping.
#'  Typically phenotype. Locus are considered independently for each group.
#' @param locus_len A \code{numeric}. The size of loci. i.e: One SNP get returned
#'  for each window of locus_len.
#'
#' @return A subset of the original \code{table} containing only the most
#' significant SNP at each locus.
#'
#' @export
topTableSNP <- function(table, pos, pval, group, locus_len = 1e6){
  tab2 <- table[, topSNP(get(pos), get(pval), locus_len = locus_len),  by = group]
  ret <- table[tab2, .SD[i.V1], on = group, by = .EACHI]
  return(ret)
}


## TODO: This only works for genic SNP and makes the package longer to load
## Find position from RSID
##
## @import SNPlocs.Hsapiens.dbSNP144.GRCh37
## @import TxDb.Hsapiens.UCSC.hg19.knownGene
## @importFrom BSgenome snpsById
## @importFrom GenomicFeatures genes
## @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle
annotateSNP <- function(rsids){
  snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  gp <- snpsById(snps, rsids, ifnotfound = "drop")
  seqlevelsStyle(gp) <- "UCSC" # "chrX"

  dt <- data.table(gp)
  return(dt)
}
