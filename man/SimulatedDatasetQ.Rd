\name{SimulatedDatasetQ}
\alias{SimulatedDatasetQ}
\docType{data}
\title{
Marker phenotypes and trait values of offspring individuals in the mapping population
}
\description{
An example data frame generated under a quadrivalent model, including marker phenotypes (i.e. SNP dosage) and trait values of mapping population for demonstration purpose. Simulation model involved a single chromosome of 100 cM long with 20 bi-allelic markers for which dosage information if known. The coefficient of double reduction at the first locus was zero. A single QTL is simulated at 34 cM (alpha=0.15 at the QTL) with parental genotype configuraion of QQqq*qqqq.
}
\usage{data("SimulatedDatasetQ")}
\format{
  A data frame with 300 observations on the following 21 variables.
  \describe{
    \item{\code{L1}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 1st locus}
    \item{\code{L2}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 2nd locus}
    \item{\code{L3}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 3rd locus}
    \item{\code{L4}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 4th locus}
    \item{\code{L5}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 5th locus}
    \item{\code{L6}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 6th locus}
    \item{\code{L7}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 7th locus}
    \item{\code{L8}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 8th locus}
    \item{\code{L9}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 9th locus}
    \item{\code{L10}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 10th locus}
    \item{\code{L11}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 11th locus}
    \item{\code{L12}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 12th locus}
    \item{\code{L13}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 13th locus}
    \item{\code{L14}}{A numeric vector indicates the allele dosage information of the 300 offspring individual on the 14th locus}
    \item{\code{L15}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 15th locus}
    \item{\code{L16}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 16th locus}
    \item{\code{L17}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 17th locus}
    \item{\code{L18}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 18th locus}
    \item{\code{L19}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 19th locus}
    \item{\code{L20}}{A numeric vector indicates the allele dosage information of the 300 offspring individuals on the 20th locus}
    \item{\code{Trait}}{A numeric vector indicates the trait value of the 300 offspring individuals}
  }
}
\details{
Here SNP alleles are denoted by integer 1 and 2 and the allele dosage indictes the copy number of allele 2 on the locus
}

\examples{
data(SimulatedDatasetQ)
}
\keyword{datasets}
