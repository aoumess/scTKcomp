\name{scTK_QnDuMAP}
\alias{scTK_QnDuMAP}
\title{Perform a "quick n' dirty" 2D uMAP on an early single cell matrix, for visualization purpose}
\usage{
  scTK_QnDuMAP <- function(in_rds = NULL, exp_name = NULL, assay = 'SLN', pca_comp = 10, my_seed = 1337, return_data = FALSE)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{exp_name}{Name of the experiment to use from the SCE}
  \item{assay}{Name of the assay in the provided experiment name, that consists in the input matrix to visualize through a uMAP reduction}
  \item{pca_comp}{Number of components to use from the intermediate PCA}
  \item{my.seed}{Seed to be used by the RNG. Systematically use a fixed value to get reproducible results.}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs a "quick n' dirty" 2D uMAP dimension reduction on an early single cell matrix from a SCE object in a RDS archive, for visualization, using default values from the Seurat methods. The modified SCE object, in which the uMAP has been added, is saved in a RDS archive with the same name as the input one ("in_rds")
}
\author{
  Bastien Job
}
