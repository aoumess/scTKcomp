\name{scTK_SeuratUMAP}
\alias{scTK_SeuratUMAP}
\title{Perform a uMAP dimension reduction using Seurat, from a reduced dimension assay of a SCE object from a RDS archive}
\usage{
  scTK_SeuratUMAP <- function(in_rds = NULL, in_dimred = NULL, ndim_max = 10, my_seed = 12345, out_rds = 'auto', return_data = FALSE, ...)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{in_dimred}{Name of the dimensionr reduction object to use from the SCE, as input for the uMAP}
  \item{ndim_max}{Number of components to use from the input dimension reduction "in_dimred"}
  \item{my.seed}{Seed to be used by the RNG. Systematically use a fixed value to get reproducible results.}
  \item{out_rds}{Path to the output RDS. If NULL, nothing will be written on disk. If 'auto', the RDS will be written along the in_rds}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs a uMAP dimension reduction using Seurat methods, from a (PCA/ICA) dimension reduction. (This function is mainly there to bypass a bug from singleCellTK that impeaches the loading of a RDS-saved SCE object containing a UMAP generated on a FMNN assay)
}
\author{
  Bastien Job
}
