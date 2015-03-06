\name{EBSeqHMM-package}
\alias{EBSeqHMM-package}
\alias{EBSeqHMM}
\docType{package}
\title{
EBSeqHMM: A Bayesian approach for identifying gene-expression changes in ordered RNA-seq experiments
}
\description{
The EBSeqHMM package implements an auto-regressive hidden Markov model for statistical analysis in ordered RNA-seq experiments (e.g. time course or spatial course data). The EBSeqHMM package provides functions to identify genes and isoforms that have non-constant expression profile over the time points/positions, and cluster them into expression paths.
}
\details{
\tabular{ll}{
Package: \tab EBSeqHMM\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2014-09-16\cr
License: \tab Artistic-2.0 \cr
}
}
\author{
Ning Leng, Christina Kendziorski
Maintainer: Ning Leng <nleng@wisc.edu>
}
\references{
Leng, N., Li, Y., Mcintosh, B. E., Nguyen, B. K., Duffin, B., Tian, S., Thomson, J. A., Colin, D., Stewart, R. M., and Kendziorski, C.
(2014). Ebseq-hmm: A bayesian approach for identifying gene-expression changes in ordered rna-seq experiments.}
\keyword{ package }
\seealso{
EBSeq
}
\examples{
data(GeneExampleData)
CondVector <- rep(paste("t",1:5,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
Sizes <- MedianNorm(GeneExampleData)
EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
           UpdateRd=2)
}
