#' @title Plot expression of a single gene
#' @usage PlotExp(NormalizedData, Conditions, Name)
#' @param NormalizedData Expression data after adjusting for library size factors
#' @param Conditions sample conditions
#' @param Name name of the gene/isoform of interest
#' @author Ning Leng
#' @examples 
#' data(GeneExampleData)
#' CondVector <- rep(paste("t",1:5,sep=""),each=3)
#' Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
#' Sizes <- MedianNorm(GeneExampleData)
#' NormData <- GetNormalizedMat(GeneExampleData, Sizes)
#' PlotExp(NormData, Conditions, "Gene_1")
#' @return PlotExp() funtion will generate line plots for genes or isoforms of interest.
#' @details PlotExp() function will generate line plots for genes or isoforms of interest.

PlotExp <- function(NormalizedData, Conditions, Name){
	GeneName <- Name
	Data <- NormalizedData
	plot(as.numeric(Conditions), Data[GeneName,],pch=20, xlab="Conditions",ylab="Expression",
     col="blue",xaxt="none",main=GeneName)
	axis(side=1, at=1:nlevels(Conditions), levels(Conditions),las=3)
	Med=sapply(levels(Conditions),function(i)median(Data[GeneName,which(Conditions==i)]))
	lines(Med, col="blue",lwd=2)

}
