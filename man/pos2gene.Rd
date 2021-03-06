% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations.R
\name{pos2gene}
\alias{pos2gene}
\title{Find genes overlapping a position}
\usage{
pos2gene(target, interval)
}
\arguments{
\item{target}{a \code{data.table} with two columns, ID and positions to match.}

\item{interval}{a \code{data.table} with at least three columns "GENE",
"start" and "end".}
}
\value{
A \code{data.table} with all positions that were within the given
 intervals along with the gene that covers them and its coordinates
}
\description{
Find genes overlapping a position
}
\note{
When there are multiple matches in the \code{interval}, this function will
return all genes overlapping the query.
}
