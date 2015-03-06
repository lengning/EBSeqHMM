\name{GeneExampleData}
\alias{GeneExampleData}
\docType{data}
\title{
	Simulated gene level data set with 5 ordered conditions
}
\description{
	'GeneExampleData' gives the gene level simulated data with 5 ordered conditions, triplicates for each condition.
	The data set was simulated following the Negative Binomial distribution. The parameters of each gene (mean and overdispersion)
	were sampled from the empirical estimates from an empirical RNA-Seq data set from Thomson lab at Morgridge Institute for 
	Research.  
}
\format{
	GeneExampleData is a matrix with 100 genes (rows) and 15 samples (columns). 
}
\seealso{
IsoExampleList
}
\examples{
data(GeneExampleData)
str(GeneExampleData)
}
\keyword{datasets}




