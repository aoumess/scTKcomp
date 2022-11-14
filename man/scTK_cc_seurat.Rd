\name{scTK_cc_seurat}
\alias{scTK_cc_seurat}
\title{Estimate the cell cycle score / phase of cells from a SCE object in a RDS archive, using Seurat}
\usage{
  scTK_cc_seurat <- function(in_rds = NULL, assay = 'counts', cc_seurat_file = NULL, nbin = 24, my_seed = 1337, out_rds = 'auto', return_data = FALSE)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{assay}{Name of the assay in the provided experiment name, that consists in the input matrix to compute the cell cycle scores / phases}
  \item{cc_seurat_file}{Path to a RDS archive that contains two genesets characterizing the S and the G2M phases}
  \item{nbin}{see ?Seurat::CellCycleScoring}
  \item{my.seed}{Seed to be used by the RNG. Systematically use a fixed value to get reproducible results.}
  \item{out_rds}{Path to the output RDS. If NULL, nothing will be written on disk. If 'auto', the RDS will be written along the in_rds}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs an estimation of the cell cycle scores / phase status of cells from a SCE object in a RDS archive, and outputs the modified SCE object containing the scores and cycle phase status as additional cell annotations in a new RDS archive.
}
\author{
  Bastien Job
}
