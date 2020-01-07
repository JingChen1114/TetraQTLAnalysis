\name{OrthogonalModel}
\alias{OrthogonalModel}
\title{
Orthogonal contrast based genetics model for autotetraploid speices
}
\description{
The orthogonal model directly dissects genetic effects into monogenic,digenic,trigenic and quadrigenic effects. This dissects genetic variance into independent components and gives consistent estimates of genetic effects provided that tetrasomic gene segregation is considered.
}
\usage{
OrthogonalModel(P1, P2, G, alpha)
}
\arguments{
  \item{P1}{
A vector(4) indicates the genotype of one parent with integer 1/2 representing trait phenotype decreasing/increasing allele.
}
  \item{P2}{
A vector(4) indicates the genotype of another parent with integer 1/2 representing trait phenotype decreasing/increasing allele.
}
  \item{G}{
A vector(5) indicates the genotypic value of genotype 1111,2111,2211,2221 and 2222 (i.e.the genotypic value is set to be zero if the corresponding genotype is not applicable). 
}
  \item{alpha}{
A numberic value indicates the coefficient of double reduction at the QTL under consideration. alpha is always equal to zero for bivalent pairing model. 
}
}

\value{
\item{Geffects}{A vector indicates population mean and genetic effects including monogenic,digenic,trigenic and quadrigenic effects if applicable.}
\item{Gvariance}{A vector indicates the component genetic variance corresponding to monogenic,digenic,trigenic and quadrigenic effects if applicable.}
}
\references{
Chen J, et al. 2018. Orthogonal contrast based models for quantitative genetic analysis in autotetraploid species. New Phytologist 220:332-346. 
}
\author{
Jing Chen
}

\examples{
P1 <- c(2,2,1,1)
P2 <- c(1,1,1,1)
alpha=0.15
G <- c(-1.035,12.971,17.979,0.0,0.0)
out <- OrthogonalModel(P1,P2,G,alpha)  
}
