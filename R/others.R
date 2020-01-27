#' dcast to a matrix
#'
#' Reformat a long data.table to a wide matrix, preserving rownames and type
#'
#' @param table The \code{data.table} to dcast.
#' @param ... Extra arguments to be passed to \code{dcast.data.table}
#'
#' @export
dcast2mat <- function(table, ...){
  # dcast the table
  dat <- dcast.data.table(data = table, ...)
  rn <- dat[[1]]
  dat <- as.matrix(dat[, -1])
  rownames(dat) <- rn
  return(dat)
  # extract rownames and column names
  # cast as matrix | this must be done after removing the rownames column to prevent additiohnal coercion
  # set the correct type
}


# @examples
# idx <- sample(1:5, 12, replace = t)
# dt1 <- data.table(original = letters[11:15][idx], expected = letters[idx])
# dt_match(dt1, letters[11:15], letters[1:5])
