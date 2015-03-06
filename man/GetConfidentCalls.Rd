\name{GetConfidentCalls}
\alias{GetConfidentCalls}
\title{Obtain confident gene calls for classifying genes into expression paths}
\usage{
GetConfidentCalls(EBSeqHMMOut, FDR=.05, cutoff=0.5, OnlyDynamic=TRUE,Paths=NULL)
}
\arguments{
  \item{EBSeqHMMOut}{output from EBSeqHMMTest function}

  \item{FDR}{Target FDR, default is 0.05.}

  \item{cutoff}{cutoff to use for defining a confident
  call. Genes with PP_path greater or equal to cutoff will
  be called as a confident call. Default is 0.5.}

  \item{OnlyDynamic}{if specifies as T, only dynamic paths
  will be shown}

  \item{Paths}{paths that are of interest. Default is NULL.
  If it is not specified, all possible paths will be
  considered.}
}
\value{
Overall: a list of genes/isoforms that are identified as DE
under the target FDR, shown are their names and PPs;
EachPath: a list object, each sublist contains confident
calls (genes/isoforms) that have PP(path)>=cutoff for a
particular expression path, shown are their names and PPs;
NumEach: length of each sublist in EachPath. EachPathName:
gene/isoform names in each of the sublists in EachPath
}
\description{
Obtain confident gene calls for classifying genes into
expression paths
}
\details{
Function GetConfidentCalls() can be used to obtain a list
of DE genes/isoforms with user specific cutoffs. To obtain
a list of DE genes/isoforms with a target FDR alpha, the
user may specify FDR=alpha. To further choose
genes/isoforms with high posterior probability of being its
most likely path, the user may specify the option cutoff
(default is 0.5). Then genes or isoforms with PP(most
likely path ) > = 0.5 will be selected
}
\note{
Output: output a list of genes that are classified to a
expression path as a confident assignment.
}
\examples{
data(GeneExampleData)
CondVector <- rep(paste("t",1:5,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
Sizes <- MedianNorm(GeneExampleData)
EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
          UpdateRd=2)
GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)
GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=TRUE)
}
\author{
Ning Leng
}

