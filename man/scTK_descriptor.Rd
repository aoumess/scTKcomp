\name{scTK_descriptor}
\alias{scTK_descriptor}
\title{Describing the content of a SingleCellExperiment object in a RDS archive}
\usage{
  scTK_descriptor <- function(in_rds = NULL, describe = 'all', sparse_level = TRUE, max_levels = 100, return_data = FALSE)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{describe}{A character specifying which kind of entry to describe. "assay" to restrict to experiments and their respective assays. "dimred" for dimension reduction objects. "coldata" for the barcode annotations. "all" for all entries}
  \item{sparse_level}{Compute the sparsity level (ie, amount of 0 values) for the 'count' assay of the main experiment.}
  \item{max_levels}{Maximum number of levels to consider a barcode annotation as a factor.}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function describes the content of a SCE object from a RDS archive.
}
\author{
  Bastien Job
}
