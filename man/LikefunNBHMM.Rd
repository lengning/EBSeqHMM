\name{LikefunNBHMM}
\alias{LikefunNBHMM}
\title{Likelihood function of the Beta-Negative Binomial HMM Model}
\usage{
LikefunNBHMM(ParamPool, InputPool)
}
\arguments{
  \item{ParamPool}{The parameters that will be estimated in
  EM.}

  \item{InputPool}{The control parameters that will not be
  estimated in EM}
}
\value{
optimal estimates of the parameters of interest
}
\description{
Likelihood function of the Beta-Negative Binomial HMM Model
}
\details{
The likelihood function of the Beta-Negative Binomial HMM
model used in EBSeqHMM. EBSeqHMM uses optim() function to
obtain the optimal estimates that minimizes the likelihood.
}
\examples{
data(GeneExampleData)
tmp <- GeneExampleData[1:10,]
In <- list(tmp,1,5,10,3,tmp,rep(1,15),as.factor(rep(1:5,each=3)), 10,cbind(rep(.5,10),rep(1,10),rep(2,10)))
Start <- c(1,1)
LikefunNBHMM(Start,In)
}
\author{
Ning Leng
}

