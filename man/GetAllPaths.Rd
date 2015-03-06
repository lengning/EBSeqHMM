\name{GetAllPaths}
\alias{GetAllPaths}
\title{Obtain all possible gene paths for an RNA-seq experiments with ordered conditions}
\usage{
GetAllPaths(EBSeqHMMOut, OnlyDynamic=TRUE)
}
\arguments{
  \item{EBSeqHMMOut}{output from EBSeqHMMTest function}

  \item{OnlyDynamic}{if specifies as TRUE, only dynamic
  paths will be shown}
}
\value{
output: a vector of paths. For example, Up-Up-Up-Up,
Up-Up-EE-EE, Up-Down-Up-EE, etc.
}
\description{
Obtain all possible gene paths for an RNA-seq experiments
with ordered conditions
}
\details{
GetAllPaths() function may be used to generate all possible
expression paths of a particular design.
}
\examples{
data(GeneExampleData)
CondVector <- rep(paste("t",1:5,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
Sizes <- MedianNorm(GeneExampleData)
EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
          UpdateRd=2)
AllPaths <- GetAllPaths(EBSeqHMMGeneOut)
}
\author{
Ning Leng
}

