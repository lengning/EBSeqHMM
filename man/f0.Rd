\name{f0}
\alias{f0}
\title{Calculate the prior predictive distribution of the Beta-Negative Binomial  model}
\usage{
f0(Input, AlphaIn, BetaIn, EmpiricalR, NumOfGroups, log)
}
\arguments{
  \item{Input}{expression values}

  \item{AlphaIn,BetaIn,EmpiricalR}{The parameters estimated
  from last iteration of EM.}

  \item{NumOfGroups}{How many transcripts within each Ng
  group}

  \item{log}{If set as TRUE, the output will in log scale.}
}
\value{
output a numeric vector, each element shows the prior
predictive probability of one gene/isoform
}
\description{
Calculate the prior predictive distribution of the
Beta-Negative Binomial model
}
\details{
Function f0() will calculate the Beta-Negative Binomial
prior predictive probability for a given set of parameters
}
\examples{
f0(matrix(rnorm(100,100,1),ncol=10), .5, .6,
           matrix(rnorm(100,200,1),ncol=10), 100, TRUE)
}
\author{
Ning Leng
}

