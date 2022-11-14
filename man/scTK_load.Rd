\name{scTK_load}
\title{Loading single cell data from Alevin/Bustools/CellRanger/UMI-tools, saving to a SCE object}
\usage{
  scTK_load <- function(data_path = NULL, sample_name = 'SAMPLE', exp_name = NULL, out_rds = 'auto', return_data = FALSE)
}
\arguments{
  \item{data_path}{Path to single cell data, be it from Alevin / Bustools / CellRanger / UMI-tools.}
  \item{sample_name}{A name to identify the current data in the ouput SCE. Please avoid non-standard characters.}
  \item{exp_name}{Name of the main experiment of the output SCE object where the assaywill be stored. An experiment name is not mandatory}
  \item{out_rds}{Path to the output RDS. If NULL, nothing will be written on disk. If 'auto', the RDS will be written along the in_rds}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs the loading into R of single cell data from multiple formats, and the outputting of the filtered barcode matrix or a disk-written SingleCellExperiment object.
}
\author{
  Bastien Job
}
