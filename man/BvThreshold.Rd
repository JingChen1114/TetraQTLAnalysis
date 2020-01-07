\name{BvThreshold}
\alias{BvThreshold}
\title{
Significance of QTLs by BvMethod
}
\description{
BvThreshold works out the threshold of LOD scores for declaring a significant QTL at 0.05 level by repeated permutation tests.
}
\usage{
BvThreshold(n, mloci, ftype, mtype, rec, o, trait, fqtl, mqtl, ir)
}
\arguments{
  \item{n}{
An integer indicates the number of offspring individuals in the segregation population.
}
  \item{mloci}{
An integer indicates the number of marker loci in the linkage group.
}
  \item{ftype}{
An integer matrix(mloci×4) indicates paternal phased genotypes on the marker loci.
}
  \item{mtype}{
An integer matrix(mloci×4) indicates maternal phased genotypes on the marker loci.
}
  \item{rec}{
A vector(mloci) indicates the recombination frequencies between the marker loci.rec[i] is the recombination frequency between the (i-1)th and ith marker and rec[1] is set to be zero.
}
  \item{o}{
An integer matrix(n×mloci×8) indicates the marker phenotypes of offspring. Each element takes a value of 1-4 representing the dosage of the corresponding allele (i.e.there are up to 8 different alleles in a outbred population) and takes a value of 0 representing the absence of the corresponding allele.
}
  \item{trait}{
A vector(n) indicates the trait phenotype data of offspring individuals in the segregation population.
}
  \item{fqtl}{
A vector(4) indicates the predicted paternal QTL configuration.
}
  \item{mqtl}{
A vector(4) indicates the predicted maternal QTL configuration.
}
  \item{ir}{
An integer indicates the number of permutation tests.
}
}

\value{

\item{threshold}{The threshold obtained for declaring a significant QTL.}

}
\references{
Chen J, et al. 2019 Mapping Quantitative Trait Loci in Autotetraploid Species.
}
\author{
Jing Chen
}
\examples{
data("LinkageMap")
data("SimulatedDatasetB")
ftype <- cbind(LinkageMap$P11,LinkageMap$P12,LinkageMap$P13,LinkageMap$P14)
mtype <- cbind(LinkageMap$P21,LinkageMap$P22,LinkageMap$P23,LinkageMap$P24)
rec <- LinkageMap$rec
mloci <- length(rec)
dosage <- as.matrix(SimulatedDatasetB[1:mloci])
colnames(dosage) <- NULL
trait <- SimulatedDatasetB$Trait
n <- length(trait)
o <- array(rep(0,n*mloci*8),dim=c(n,mloci,8))
for(i in 1:n){                 
   for(j in 1:mloci){
      o[i,j,1]=dosage[i,j]
      o[i,j,2]=4-dosage[i,j]
   }
}
fqtl <- c(2,2,1,1)
mqtl <- c(1,1,2,1)
ir <- 100
out <- BvThreshold(n, mloci, ftype, mtype, rec, o, trait, fqtl, mqtl, ir)
}