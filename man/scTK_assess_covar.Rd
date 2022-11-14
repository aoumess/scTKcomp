\name{scTK_assess_covar}
\alias{scTK_assess_covar}
\title{Perform a quick estimation of the "weight" of cell annotations on a single cell data matrix}
\usage{
  scTK_assess_covar <- function(in_rds = NULL, exp_name = NULL, assay = 'counts', factor_names = NULL, conti_names = NULL, ctrl_features = NULL, marker_features = NULL, red_method = 'pca', ndim_max = 10, center = TRUE, scale = TRUE, coef_cut = 0, out_png = 'auto', color_palette = c("white", "orangered3"))
}
\arguments{
  \item{in_rds}{Path to a RDS archive that contains a SCE object}
  \item{exp_name}{Name of the experiment to use from the SCE}
  \item{assay}{Name of the assay containing the matrix to use from the provided experiment name}
  \item{factor_names}{Name of one (or more) factor cell annotations from the SCE object to assess}
  \item{conti_names}{Name of one (or more) continuous cell annotations from the SCE object to assess}
  \item{ctrl_features}{Name of one (or more) features (ie, genes) to use as a control (which expression should have no correlation with the assay)}
  \item{marker_features}{Name of one (or more) features (ie, genes) to use as a markers (which expression is expected to have correlation with the assay)}
  \item{red_method}{Dimension reduction method to use for the assessment. Should be 'pca', 'mds.euc', or 'mds.pear'}
  \item{ndim_max}{Number of components to use from the intermediate dimension reduction}
  \item{center}{Centers the data before the dimension reduction}
  \item{scale}{Scales the data before the dimension reduction}
  \item{coef_cut}{Do not display any color in the resulting heatmap when a score is lower than this value}
  \item{out_png}{Path to the output PNG image containing the heatmap. If 'auto', the PNG will be written along the "in_rds" archive}
  \item{color_palette}{Vector of twho color names to generate the color palette of the heatmap. The first and second values correspond to lower and higher scores, respectively}
}
\description{
  This function attempts to assess the "weight" (by a correlation test coefficient for contiunous variables, a Kruskal-Wallis test statistic for factors) of covariates in a single cell expression matrix. The output is a heatmap showing weights as coloured cells.
}
\author{
  Bastien Job
}
