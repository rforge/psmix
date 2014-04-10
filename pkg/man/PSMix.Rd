\name{PSMix}
\alias{PSMix}
\title{ Population structure inference using likelihood-based approach }
\description{
  Population structure inference using likelihood-based approach
}
\usage{
  PSMix(K=2, Geno, itMax=10000, eps=1e-06, seed=NULL, MarkerVar=FALSE, verbose=FALSE)
}
\arguments{
  \item{K}{ number of populations }
  \item{Geno}{ the genotype data }
  \item{itMax}{ maximum number of iterations }
  \item{eps}{ convergence precision }
  \item{seed}{ random seed }
  \item{MarkerVar}{ whether marker dependence is allowed }
  \item{verbose}{ print out the maximum likelihood information }
}
\value{
  \item{AmPr}{estimated admixture proportions for each sample}
  \item{AmId}{estimated underlying population membership based on the admixture
proportions for each sample}
}
\references{ 
  Tang, H. and Peng, J. and Wang, P. and Risch, N. (2005)
   \emph{Estimation of individual admixture: analytical and study design 
considerations.}
   \bold{Genetic Epidemiology}, 28, 289-301.
  
  Wu, B. and Liu, N. and Zhao, H. (2006)
   \emph{PSMIX: an R package for population structure inference via maximum 
likelihood method.}
   \bold{BMC Bioinformatics}, 7:317.
}
\seealso{
  see \code{\link{geno.ps.x}} for the genotype data format
}
\examples{
\dontrun{
 ## Load in the Pima-Surui Data"
  data(geno.ps.x)
  data(geno.ps.y)
 ## Infer population structure with first 50 markers
 obj1 = PSMix(K=2, geno.ps.x[,1:50])
 ## admixture probability estimation
 round(obj1$AmPr, 3)
 ### Plot the probability estimation: Pima/Surui are coded using different symbols/colors.
 ##   And samples are sorted according to the probabilities.
 i1 = order(obj1$AmPr[,1])
 t1 = as.integer(as.factor(geno.ps.y[i1]))
 plot(obj1$AmPr[i1,1], pch=t1,col=t1, xlab="Samples", ylab="Prob" )
 abline(h=0.5, col=3, lty=2)
 ## admixture membership estimation
 table(geno.ps.y, obj1$AmId)
 ## Infer population structure with first 100 markers
 obj2 = PSMix(K=2, geno.ps.x[,1:100])
 ## admixture probability estimation
 round(obj2$AmPr,3)
 i2 = order(obj2$AmPr[,1])
 t2 = as.integer(as.factor(geno.ps.y[i2]))
 ### Plot the probability estimation: Pima/Surui are coded using different symbols/colors.
 ##   And samples are sorted according to the probabilities.
 dev.new()
 plot(obj2$AmPr[i2,1], pch=t2,col=t2, xlab="Samples", ylab="Prob" )
 abline(h=0.5, col=3, lty=2)
 ## admixture membership estimation
 table(geno.ps.y, obj2$AmId)
}
}
\keyword{ function }
