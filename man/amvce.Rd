\name{amvce}
\alias{amvce}
\title{
Linear mixded model (or quantitative genetic animal model) variance component estimation using REML
}
\description{
Linear mixded model (or quantitative genetic animal model) variance component estimation using REML. Author: Minghui Wang <m.h.wang@live.com>
}
\usage{
amvce(y, X, A, tolerance = 0.001, maxIter = 100, verbose = FALSE)
}
\arguments{
  \item{y}{
a vector of numeric values
}
  \item{X}{
a matrix of fixed effects
}
  \item{A}{
a square matrix of random effect correlation structure
}
  \item{tolerance}{
(optional) a relative tolerance to test converging
}
  \item{maxIter}{
(optional) maximum number of iteration
}
  \item{verbose}{
(optional) whether to print progress
}
}

\value{
\item{coefficients}{fixed effect coefficients}
\item{VC}{a vector of variance components}
\item{loglik}{log-likelihood}
}
\references{
Yang S, Liu Y, Jiang N, Chen J, Leach L, Luo Z and Wang M. Genome-wide eQTLs and heritability for gene expression traits in unrelated individuals. BMC Genomics 2014, 15:13
}
\examples{

}
