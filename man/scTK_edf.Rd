\name{scTK_edf}
\alias{scTK_edf}
\title{Filtering empty droplets from a SCE object in a RDS archive, using DropletUtils::emptyDrops}
\usage{
  scTK_edf <- function(in_rds = NULL, assay = 'counts', droplets_min = 1E+04, emptydrops_fdr = 1E-03, emptydrops_retain = NULL, my_seed = 1337, draw_plots = TRUE, out_rds = 'auto', return_data = FALSE)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{assay}{Name of the assay in the main experiment, that consists in the input matrix for DropletUtils::emptyDrops}
  \item{droplets_min}{Minimum number of expected droplets. If the actual content in the input assay is lower than this value, no filtering is performed. This is to avoid inadvertedly trying a second pass by mistake}
  \item{emptydrops_fdr}{Maximum FDR-adjusted p-value level from the emptyDrops result to consider a droplet as non-empty. See \code{'?DropletUtils::emptyDrops'}, section "Value"}
  \item{emptydrops_retain}{Number of cells. Value to force the retention of a fixed number of cells. The NULL default value means that the value used will be the cut-off determined by the \code{emptyDrops::barcodeRanks} function.}
  \item{my.seed}{Seed to be used by the RNG. Systematically use a fixed value to get reprocible results.}
  \item{draw_plots}{Outputs the famous kneeplot.}
  \item{out_rds}{Path to the output RDS. If NULL, nothing will be written on disk. If 'auto', the RDS will be written along the in_rds}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs the filtering of empty droplets using the emptyDrops method from \code{DropletUtils}, and the output of modified SCE object to a variable or ina disk-written RDS archive.
}
\author{
  Bastien Job
}
