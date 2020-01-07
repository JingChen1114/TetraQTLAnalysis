\name{TraitSeg}
\alias{TraitSeg}
\title{
Detection of major gene segregation in an outbred autotetraploid population
}
\description{
The segregation analysis models the trait phenotype data distribution as a mixed normal distribution in which each component distribution corresponds to a particular genotype of the major QTL. THe trait segregation analysis is carried out across all twelve possible parental genotype configuraions, across the range of possible value (0-0.25) for the rate of double reduction. An EM algorithm is implemented to calculate the MLEs of model parameters and a likelihood ratio test is applied to test the significance of major QTL segregation.
}
\usage{
TraitSeg(n, trait)
}
\arguments{
  \item{n}{
An integer indicates the number of individuals in the segregation population.
}
  \item{trait}{
A vector(n) indicates the trait phenotype data of individuals in the segregation population.
}
}

\value{
\item{p1}{Estimated parental genotype of parent 1}
\item{p2}{Estimated parental genotype of parent 2}
\item{alpha}{MLE of the coefficient of double reduction}
\item{g_value}{A vector(5) indicates MLEs of the genotypic value for genotype qqqq,Qqqq,QQqq,QQQq and QQQQ.}
\item{p_value}{P value of the significance of the major QTL segregation}
\item{variance}{Total variance of the trait phenotype data}
\item{g_var}{Genetic variance contributed by the major gene}
}
\references{
Chen J, et al. 2018. Orthogonal contrast based models for quantitative genetic analysis in autotetraploid species. New Phytologist 220:332-346. 
}
\author{
Jing Chen
}

\examples{
data("TraitData")
trait <- as.matrix(TraitData)
n <- nrow(trait)
ir <- ncol(trait)/2
trait1 <- rep(0,n)
for(i in 1:n){
  s=0
  ir1=ir
  for(j in 1:ir){
    if(is.na(trait[i,j*2])){
      ir1=ir1-1
    }else{
      s=s+trait[i,j*2]
    }
  }
  trait1[i]=s/ir1
}
out <- TraitSeg(n,trait1)    
}

