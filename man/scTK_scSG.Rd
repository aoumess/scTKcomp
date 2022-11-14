\name{scTK_scSG}
\alias{scTK_scSG}
\title{Identify "sensitive genes", rather than usual "highly variable genes" for a reduction dimension of a SCE object by selection at the feature level}
\usage{
  scTK_scSG <- function(in_rds = NULL, exp_name = NULL, raw_assay = 'counts', norm_assay = 'SLN', n_features = 2000, n_PC = 10, resolution = 0.6, out_rds = 'auto', return_data = FALSE)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{exp_name}{Name of the experiment to use from the SCE}
  \item{raw_assay}{Name of the assay containing the raw count matrix in the provided experiment name}
  \item{norm_assay}{Name of the assay containing the normalized matrix in the provided experiment name}
  \item{n_feature}{Number of "sensitive genes" to retrieve}
  \item{n_PC}{Number of components to use from the intermediate PCA}
  \item{resolution}{Resolution to use for the Louvain clustering. This consists in a rough clustering, thus a low value is recommended.}
  \item{out_rds}{Path to the output RDS. If NULL, nothing will be written on disk. If 'auto', the RDS will be written along the in_rds}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs the identification of the "top n_feature" sensitive genes, using scSensitiveGenes, from a SCE object in a RDS archive.
}
\author{
  Bastien Job
}
