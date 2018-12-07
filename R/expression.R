#' Process kallisto outputs
#'
#' @param bstrapdir A \code{character}. The path of a kallisto bootstrap directory.
#' @param quantilenorm A \code{logical}. If set to TRUE, the data is quantile
#'  normalized to standardize expression across samples.
#'
#' @return
#' A \code{data.table} in long format.
#'
#'
#' @export
readKali <- function(bstrapdir, quantilenorm = F){
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
  tpms <- rbindlist(ablist)
  tpms[, transcript := gsub("\\..*", "", target_id)]
  if(quantilenorm){
    print("Quantile normalization")
    tpms <- .quantnorm_kali_dt(tpms, "transcript", "sample", "tpm")
  }
  return(tpms)
}


#' @importFrom preprocessCore normalize.quantiles
.quantnorm_kali_dt <- function(dt, v1, v2, valuevar){
  dc <- dcast(dt, get(v1) ~ get(v2), value.var=valuevar)
  mat <- as.matrix(dc[, 2:ncol(dc)])
  matqn <- normalize.quantiles(mat)
  rownames(matqn) <- dc$v1
  colnames(matqn) <- colnames(mat)
  qnormdt <- data.table(matqn, keep.rownames = T)
  qnormdt <- melt(qnormdt, id.vars="rn", variable.name = v2, value.name=valuevar)
  setnames(qnormdt, "rn", v1)
  return(qnormdt)
}

#' topTable to data.table
#'
#' Transform limma topTable output into a data.table object
#'
#' @param tt
#'
#' @return A \code{data.table} with the gene names as the first column
#'
#' @export
tt2dt <- function(tt){
  gns <- rownames(tt)
  dt <- data.table(tt)
  dt[, gene := gns][]
  setcolorder(dt, c(ncol(dt), 1:(ncol(dt)-1)))
  return(dt)
}


#' @export
txi2de <- function(txi, pheno){
  if(any(colnames(txi$counts) != rownames(pheno)))
    stop("The column names of `txi` should match the rownames of `pheno`") 
  design <- model.matrix(~gender, data = pheno) #Ensure ordering is identical (messes it up silently otherwise)
  y <- DGEList(txi$counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, ]
  y <- calcNormFactors(y) #TMM normalization
  v <- voom(y, design) #Transform counts into log2CPM
  # Fit
  fit <- lmFit(v, design = design)
  eb <- eBayes(fit)
  tt <- topTable(eb, coef="gendermale", number = Inf)
  dt <- tt2dt(tt)
}
