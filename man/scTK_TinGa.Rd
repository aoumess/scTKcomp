\name{scTK_TinGa}
\alias{scTK_TinGa}
\title{Perform a trajectory analysis of a SCE object from a RDS archive, using TinGa}
\usage{
  scTK_scSG <- function(in_rds = NULL, in_dimred = 'counts', clusters = NULL, out_dir = 'auto', return_data = FALSE)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{in_dimred}{Name of the dimensionr reduction object to use from the SCE, as input for the uMAP}
  \item{clusters}{Name of a cell annotation of the SCE object to use as cluster results}
  \item{out_dir}{Path to the output dir where the output RDS and plots will be written. If 'auto', the output files will be written in a subfolder named from in_rds}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs trajectory analysis using TinGa, from a SCE object in a RDS archive.
}
\author{
  Bastien Job
}
