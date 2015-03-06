\name{IsoExampleList}
\alias{IsoExampleList}
\docType{data}
\title{
	Simulated isoform level data set with 5 ordered conditions
}
\description{
	'IsoExampleList' gives the isoform level simulated data with 5 ordered conditions, triplicates for each condition.
  The data set was simulated following the Negative Binomial distribution. The parameters of each isoform (mean and overdispersion)
	  were sampled from the isoform level empirical estimates from an empirical RNA-Seq data set from Thomson lab at Morgridge Institute for 
		  Research.   
}
\format{
	  IsoExampleList is a list with three components. IsoExampleList$IsoExampleData contains a
		matrix with 200 isoform (rows) and 15 samples (columns). IsoExampleList$IsoNames contains a vector of isoform names.
		IsoformExampleList$IsosGeneNames contains a vector indicating the gene each isoform belongs to.
}
\seealso{
GeneExampleData
}
\examples{
data(IsoExampleList)
str(IsoExampleList)
}
\keyword{datasets}




