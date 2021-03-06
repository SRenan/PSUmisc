% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/others.R
\name{dcast2mat}
\alias{dcast2mat}
\title{dcast to a matrix}
\usage{
dcast2mat(table, ...)
}
\arguments{
\item{table}{The \code{data.table} to dcast.}

\item{...}{Extra arguments to be passed to \code{dcast.data.table}}
}
\description{
Reformat a long data.table to a wide matrix, preserving rownames and type
}
