\name{Heritability}
\alias{Heritability}
\title{
Narrow-sense heritability estimation
}
\description{
A linear mixed model is developed to estimate the narrow-sense heritability of a quantitative trait in a full-sib autotetraploid population.
}
\usage{
Heritability(ir, trait, alpha)
}
\arguments{
  \item{ir}{
An integer indicates the number of experimental replicates.
}
  \item{trait}{
A numerical matrix(n*(2ir)) indicates the trait phenotype data of a segregation population with ir replicates. The odd column indicate the individual ID and the even colunm indicate the corresponding phenotype value. Missing data is denoted by NA.
}
  \item{alpha}{
A numerical value indicates the coefficient of double reduction.
}
}
\value{
\item{V_g}{Estimated genetic variance}
\item{V_e}{Residual variance}
\item{h2n}{Estimated narrow-sense heritability}
}
\references{
Chen J, et al. 2018. Orthogonal contrast based models for quantitative genetic analysis in autotetraploid species. New Phytologist 220:332-346. 
}
\author{
Jing Chen
}
\examples{
alpha=0.055
data("TraitData")
trait <- as.matrix(TraitData)
ir <- ncol(trait)/2
out <- Heritability(ir,trait,alpha)
}
