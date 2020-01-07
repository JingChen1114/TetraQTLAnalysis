\name{TraitData}
\alias{TraitData}
\docType{data}
\title{
Trait phenotype data of a segregation population
}
\description{
An example data frame of trait phenotype data of a segregation population with three replicates. The 1st, 3rd and 5th column indicate the individual ID and the 2nd, 4th and 6th colunm indicate the corresponding phenotype value. Missing data is denoted by NA.
}
\usage{data("TraitData")}
\format{
  A data frame with 291 observations on the following 6 variables.
  \describe{
    \item{\code{N1}}{A numeric vector indicates the individual ID in the first replicate}
    \item{\code{R1}}{a numeric vector indicates the phenotype value in the first replicate}
    \item{\code{N2}}{a numeric vector indicates the individual ID in the second replicate}
    \item{\code{R2}}{a numeric vector indicates the phenotype value in the second replicate}
    \item{\code{N3}}{a numeric vector indicates the individual ID in the third replicate}
    \item{\code{R3}}{a numeric vector indicates the phenotype value in the third replicate}
  }
}

\examples{
data(TraitData)
}
\keyword{datasets}
