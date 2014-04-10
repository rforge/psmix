\name{geno.ps.x}
\alias{geno.ps.x}
\alias{geno.ps.y}
\title{ Pima Surui genotype data }
\description{
  \code{geno.ps.x} 377 autosomal microsatellites information for 25 Pima and 
21 Surui samples. two consecutive rows are used to describe the two alleles 
for one sample, each column corresponds to one microsatellite. Missing data 
is coded as NA.

  \code{geno.ps.x} sample status vector, being Pima or Surui.
}
\usage{data(geno.ps.x)}
\details{
  The Pima Suri data was reported at Rosenberg et.al 2001 science paper.
}
\source{
  \url{http://rosenberglab.bioinformatics.med.umich.edu/diversity.html}
}
\references{
  Rosenberg, Noah A. and Pritchard, Jonathan K. and Weber, James L. and Cann, 
  Howard M. and Kidd, Kenneth K. and Zhivotovsky, Lev A. and Feldman, Marcus W. (2002).
  \emph{Genetic Structure of Human Populations.}
  \bold{Science}, 298, 2381-2385.
}
\examples{
\dontrun{
  data(geno.ps.x)
  data(geno.ps.y)
  ## sample status sumary
  table(geno.ps.y)
  ## marker information
  table(geno.ps.x[,3], exclude=NULL)
  ## identify those subset of markers without missing information
  ps.xc = as.data.frame( t(na.exclude(t(geno.ps.x))) )
  ## only 113 markers left
  dim(ps.xc)
}
}
\keyword{datasets}
