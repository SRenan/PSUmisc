#' Keep only the latest log for each job
#'
#' This function finds all PBS log files (e/o) in a directory and keep only the
#' latest run for each job.
#'
#' @param logdir A \code{character}. The directory with log files to cleanup.
#'
#' @return A \code{data.table}. Summary of the remaining logs.
#'
#' @importFrom tools file_path_sans_ext
#' @export
keepLatestLog <- function(logdir){
  logs <- list.files(logdir, full.names = T, include.dirs = F)
  logdt <- data.table(file = logs, job = file_path_sans_ext(basename(logs)))
  logdt[grep("\\.e", logs), type := "e"]
  logdt[grep("\\.o", logs), type := "o"]
  logdt[, jobid := as.numeric(gsub("o|e", "", file_ext(file)))]
  # Keep maximum jobid
  # TODO: Add a check and/or option to remove by file age
  keeplog <- logdt[jobid %in% logdt[logdt[, .I[which.max(jobid)], by = c("job")]$V1]$jobid]
  rmlog  <- logdt[!jobid %in% logdt[logdt[, .I[which.max(jobid)], by = c("job")]$V1]$jobid]
  message(paste("Removing", nrow(rmlog), "files in", logdir))
  unlink(rmlog$file)
  return(keeplog)
}
