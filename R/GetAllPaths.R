#' @title Obtain all possible gene paths for an RNA-seq experiments with ordered conditions
#' @usage GetAllPaths(EBSeqHMMOut, OnlyDynamic=TRUE)
#' @param EBSeqHMMOut output from EBSeqHMMTest function
#' @param OnlyDynamic if specifies as TRUE, only dynamic paths will be shown
#' @details GetAllPaths() function may be used to generate all possible expression paths of a particular design.
#' @return output: a vector of paths. For example, Up-Up-Up-Up, Up-Up-EE-EE, Up-Down-Up-EE, etc.
#' @author Ning Leng
#' @examples 
#' data(GeneExampleData)
#' CondVector <- rep(paste("t",1:5,sep=""),each=3)
#' Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
#' Sizes <- MedianNorm(GeneExampleData)
#' EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
#'           UpdateRd=2)
#' AllPaths <- GetAllPaths(EBSeqHMMGeneOut)

GetAllPaths <- function(EBSeqHMMOut, OnlyDynamic=TRUE){

AllPaths <- colnames(EBSeqHMMOut$MgAllPP)

WithEEPath <- grep("EE",AllPaths)
DynPath <- AllPaths[-WithEEPath]

if(OnlyDynamic==TRUE)PathsConsider <- DynPath
if(OnlyDynamic==FALSE)PathsConsider <- AllPaths

PathsConsider
}

