\name{EBTest_ext}
\alias{EBTest_ext}
\title{Extented EBTest function}
\usage{
EBTest_ext(Data,NgVector=NULL,Conditions,
	sizeFactors, maxround, Pool=FALSE, NumBin=1000,
	ApproxVal=10^-10, Alpha=NULL, Beta=NULL,
	PInput=NULL,RInput=NULL,PoolLower=.25,
	PoolUpper=.75,OnlyCalcR=FALSE,Print=TRUE)
}
\arguments{
  \item{Data}{Input data, rows are genes/isoforms and
  columns are samples. Data should come from a two
  condition experiment}

  \item{NgVector}{Ng vector; NULL for gene level data}

  \item{Conditions}{A factor indicates the condition
  (time/spatial point) which each sample belongs to. Only
  two levels are allowed.}

  \item{sizeFactors}{a vector indicates library size
  factors}

  \item{maxround}{number of iteration}

  \item{Pool}{While working without replicates, user could
  define the Pool = TRUE in the EBTest function to enable
  pooling.}

  \item{NumBin}{By defining NumBin = 1000, EBSeq will group
  the genes with similar means together into 1,000 bins.}

  \item{PoolLower,PoolUpper}{With the assumption that only
  subset of the genes are DE in the data set, we take genes
  whose FC are in the PoolLower - PoolUpper quantile of the
  FCs as the candidate genes (default is 25%-75%). For each
  bin, the bin-wise variance estimation is defined as the
  median of the cross condition variance estimations of the
  candidate genes within that bin.  We use the cross
  condition variance estimations for the candidate genes
  and the bin-wise variance estimations of the host bin for
  the non-candidate genes.}

  \item{ApproxVal}{The variances of the transcripts with
  mean < var will be approximated as mean/(1-ApproxVal).}

  \item{Alpha,Beta,PInput,RInput}{If the parameters are
  known and the user doesn't want to estimate them from the
  data, user may specify them here.}

  \item{Print}{Whether print the elapsed-time while running
  the test.}

  \item{OnlyCalcR}{if OnlyCalcR=TRUE, the function will
  only return estimation of r's.}
}
\value{
See \code{\link{EBTest}}
}
\description{
Extented EBTest function
}
\details{
EBSeq_ext() function is an extension of EBTest() function,
which is used to calculate the conditional probability
P(X_g,t | X_g,t-1). In EBSeqHMM, we assume the conditional
distribution is Beta-Negative Binomial.
}
\examples{
data(GeneExampleData)
Data=GeneExampleData[,1:6]
CondVector <- rep(paste("t",1:2,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2"))
Sizes <- MedianNorm(Data[1:10,])
Out <- EBTest_ext(Data=Data[1:10,], sizeFactors=Sizes, Conditions=Conditions,
         maxround=1)
}
\author{
Ning Leng
}

