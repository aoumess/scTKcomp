\name{scTK_recompress}
\alias{scTK_recompress}
\title{Recompressing a RDS object}
\usage{
  scTK_recompress <- function(in_rds = NULL, compress = 'bzip2')
}
\arguments{
  \item{in_rds}{Path to a RDS archive to recompress.}
  \item{compress}{Compression algorithm to use. Should be "gzip", "bzip2" or "xz"}
}
\description{
  This function recompresses a RDS archive to another RDS archive, with a different compression method. "bzip2" is highly recommended for its compression ratio versus compression time.
}
\author{
  Bastien Job
}
