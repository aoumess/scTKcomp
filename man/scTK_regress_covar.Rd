\name{scTK_regress_covar}
\alias{scTK_regress_covar}
\title{Perform a regression of on or multiple factor or continuous covariates on a single cell expression matrix}
\usage{
  scTK_regress_covar <- function(in_rds = NULL, exp_name = NULL, assay = 'counts', vars_to_regress = NULL, model_use = 'linear', scale_residuals = TRUE, scale_limit = 10, center_residuals = TRUE, out_rds = 'auto', suffix = NULL, return_data = FALSE)
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{exp_name}{Name of the experiment to use from the SCE}
  \item{assay}{Name of the assay containing the matrix to use from the provided experiment name}
  \item{vars_to_regress}{Name of one (or more) factor or continuous cell annotations from the SCE object to regress}
  \item{model_use}{Type of the regression model. Should be "poisson" for a counts matrix, "linear" for other types of matrix}
  \item{scale_residuals}{Scales the residues obtained after the regression}
  \item{scale_limit}{Trim values so that their absolute is equal to scale_limit. This is to downscale outliers remaining despite scaling.}
  \item{center_residuals}{Centers the residues obtained after the regression}
  \item{out_rds}{Path to the output RDS. If NULL, nothing will be written on disk. If 'auto', the RDS will be written along the in_rds}
  \item{return_data}{Return the SCE object to a variable}
}
\description{
  This function performs the regression of one or multiple, factor or continuous covariates on a single cell matrix. The main purpose is to remove putative sources of bias.
}
\author{
  Bastien Job
}
