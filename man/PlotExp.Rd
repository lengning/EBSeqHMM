\name{PlotExp}
\alias{PlotExp}
\title{Plot expression of a single gene}
\usage{
PlotExp(NormalizedData, Conditions, Name)
}
\arguments{
  \item{NormalizedData}{Expression data after adjusting for
  library size factors}

  \item{Conditions}{sample conditions}

  \item{Name}{name of the gene/isoform of interest}
}
\value{
PlotExp() funtion will generate line plots for genes or
isoforms of interest.
}
\description{
Plot expression of a single gene
}
\details{
PlotExp() function will generate line plots for genes or
isoforms of interest.
}
\examples{
data(GeneExampleData)
CondVector <- rep(paste("t",1:5,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
Sizes <- MedianNorm(GeneExampleData)
NormData <- GetNormalizedMat(GeneExampleData, Sizes)
PlotExp(NormData, Conditions, "Gene_1")
}
\author{
Ning Leng
}

