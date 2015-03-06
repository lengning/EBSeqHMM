\name{beta.mom}
\alias{beta.mom}
\title{Method of moments estimation ( beta distribution )}
\usage{
beta.mom(qs.in)
}
\arguments{
  \item{qs.in}{A vector contains the numbers that will be
  fitted with a beta distribution.}
}
\value{
alpha.hat,beta.hat: Returns the estimation of alpha and
beta.
}
\description{
Method of moments estimation ( beta distribution )
}
\details{
beta.mom() function can be used to estimate parameters in a
Beta function using method of moments
}
\examples{
beta.mom(rbeta(10,1,1))
}
\author{
Ning Leng
}

