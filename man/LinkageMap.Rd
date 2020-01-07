\name{LinkageMap}
\alias{LinkageMap}
\docType{data}
\title{
Linkage Map 
}
\description{
A data frame including information of linkage map and parental genotypes. Simulation model involved a single chromosome of 100 cM long with 20 bi-allelic markers for which dosage information if known. The coefficient of double reduction at the first locus was zero. A single QTL is simulated at 34 cM (alpha=0.15 at the QTL).
}
\usage{data("LinkageMap")}
\format{
  A data frame with 20 observations on the following 10 variables.
  \describe{
    \item{\code{P11}}{A numeric vector indicates the paternal marker allele on the 1st chromosome on the 20 loci}
    \item{\code{P12}}{A numeric vector indicates the paternal marker allele on the 2nd chromosome on the 20 loci}
    \item{\code{P13}}{A numeric vector indicates the paternal marker allele on the 3rd chromosome on the 20 loci}
    \item{\code{P14}}{A numeric vector indicates the paternal marker allele on the 4th chromosome on the 20 loci}
    \item{\code{P21}}{A numeric vector indicates the maternal marker allele on the 1st chromosome on the 20 loci}
    \item{\code{P22}}{A numeric vector indicates the maternal marker allele on the 2nd chromosome on the 20 loci}
    \item{\code{P23}}{A numeric vector indicates the maternal marker allele on the 3rd chromosome on the 20 loci}
    \item{\code{P24}}{A numeric vector indicates the maternal marker allele on the 4th chromosome on the 20 loci}
    \item{\code{alpha}}{A numeric vector indicates the coefficient of double reduction on the 20 loci}
    \item{\code{rec}}{A numeric vector indicates the recombination frequencies between markers.rec[i] is the recombination frequency between the (i-1)th and ith marker and rec[1] is set to be zero.}
  }
}

\examples{
data(LinkageMap)
}
\keyword{datasets}
