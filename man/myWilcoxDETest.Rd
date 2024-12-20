% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forked_seurat.R
\name{myWilcoxDETest}
\alias{myWilcoxDETest}
\title{Modified Seurat::WilcoxDETest by extending the output with statistics in addition to p-values}
\usage{
myWilcoxDETest(data.use, cells.1, cells.2, verbose = TRUE, ...)
}
\arguments{
\item{data.use}{Data matrix to test}

\item{cells.1}{Group 1 cells}

\item{cells.2}{Group 2 cells}

\item{verbose}{Print a progress bar}

\item{...}{Extra parameters passed to wilcox.test}
}
\value{
Returns a p-value ranked matrix of putative differentially expressed features
}
\description{
Modified Seurat::WilcoxDETest by extending the output with statistics in addition to p-values
}
