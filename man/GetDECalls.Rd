\name{GetDECalls}
\alias{GetDECalls}
\title{Obtain DE gene/isoform list at a certain FDR}
\usage{
GetDECalls(EBSeqHMMOut,FDR=.05)
}
\arguments{
  \item{EBSeqHMMOut}{output from EBSeqHMMTest function}

  \item{FDR}{Target FDR; default is 0.05}
}
\value{
a list of genes/isoforms that are identified as DE under
the target FDR, shown are their names and PPs;
}
\description{
Obtain DE gene/isoform list at a certain FDR
}
\details{
Function GetDECalls() can be used to obtain a list of DE
genes/isoforms with user specific cutoffs. To obtain a list
of DE genes/isoforms with a target FDR alpha, the user may
specify FDR=alpha.
}
\note{
Output: output a list of genes that are DE in at least one
condition in an RNA-seq experiment with multiple ordered
conditions
}
\examples{
data(GeneExampleData)
CondVector <- rep(paste("t",1:5,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
Sizes <- MedianNorm(GeneExampleData)
EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
          UpdateRd=2)
GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)
}
\author{
Ning Leng
}

