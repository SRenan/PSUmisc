#' Transform lm/glm output to data.table
#'
#' For a list of positions, find the closest match in a given reference.
#'
#' @param pos A \code{lm} object. The output of \code{lm} or \code{glm}
#' @export
fit2dt <- function(fit){
  if(inherits(fit, "lm")){
    fit <- summary(fit)
  }
  coefs <- fit$coefficients
  form <- format(eval(fit$call[[2]]))
  if(length(coefs)==0){
    dt <- data.table()
  } else{
    dt <- data.table(coefs, keep.rownames = T)
    dt[, model := form]
    setnames(dt, c("covariate", "estimate", "serr", "tval", "pval", "model"))
  }
  return(dt)
}
