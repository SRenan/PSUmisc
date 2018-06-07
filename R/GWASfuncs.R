library(gtools)
library(ggplot2)
library(data.table)

#' Draw manhattan plot
#'
#' @param chr Vector of chromosome names
#' @param pos Vector of genomic positions
#' @param pval Vector of P-values from association test
#'
#' @importFrom ggplot2 ggplot theme theme_bw element_blank
#' @importFrom ggplot2 geom_point geom_hline geom_rect
#' @importFrom ggplot2 scale_fill_manual scale_colour_manual scale_x_continuous
#' @importFrom ggplot2 xlab ylab
#' @importFrom gtools mixedorder
#' @export
ggman <- function(chr, pos, pval, names = NULL, minpval = 1e-50, signif = 5e-8, lowp = 5e-3){
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
