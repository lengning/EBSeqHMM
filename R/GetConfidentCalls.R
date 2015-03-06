#' @title Obtain confident gene calls for classifying genes into expression paths
#' @usage GetConfidentCalls(EBSeqHMMOut, FDR=.05, cutoff=0.5, OnlyDynamic=TRUE,Paths=NULL)
#' @param EBSeqHMMOut output from EBSeqHMMTest function
#' @param FDR Target FDR, default is 0.05.
#' @param cutoff cutoff to use for defining a confident call. Genes with PP_path greater or 
#' equal to cutoff will be called as a confident call. Default is 0.5.
#' @param OnlyDynamic if specifies as T, only dynamic paths will be shown
#' @param Paths paths that are of interest. Default is NULL. If it is not specified, all possible 
#' paths will be considered.
#' @note Output: output a list of genes that are classified to a expression path as a confident assignment.
#' @author Ning Leng
#' @examples 
#' data(GeneExampleData)
#' CondVector <- rep(paste("t",1:5,sep=""),each=3)
#' Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
#' Sizes <- MedianNorm(GeneExampleData)
#' EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
#'           UpdateRd=2)
#' GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)
#' GeneConfCalls <- GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=TRUE)
#' @details Function GetConfidentCalls() can be used to obtain a list of DE genes/isoforms
#' with user specific cutoffs. To obtain a list of DE genes/isoforms with
#' a target FDR alpha, the user may specify FDR=alpha. To further choose
#' genes/isoforms with high posterior probability of being its most likely path,
#' the user may specify the option cutoff (default is 0.5). Then genes or isoforms with PP(most likely path ) > = 0.5 will be selected
#' @return Overall: a list of genes/isoforms that are identified as DE under the target FDR, shown are their names and PPs;
#' EachPath: a list object, each sublist contains confident calls 
#' (genes/isoforms) that have PP(path)>=cutoff for a particular expression path, shown are their names and PPs;
#' NumEach: length of each sublist in EachPath.
#' EachPathName: gene/isoform names in each of the sublists in EachPath

GetConfidentCalls <- function(EBSeqHMMOut,FDR=0.05, cutoff=.5, OnlyDynamic=TRUE,Paths=NULL){

SigPPLarge <- EBSeqHMMOut$MgAllPP
SigMAPLargeChar <- EBSeqHMMOut$MgAllMAPChar
SigMaxValLarge <- EBSeqHMMOut$MgAllMaxVal

AllPaths <- colnames(EBSeqHMMOut$MgAllPP)
WithUp <- grep("Up",AllPaths)
WithDown <- grep("Down",AllPaths)
UpAndDown <- union(WithUp, WithDown)
AllEEPath <- AllPaths[-UpAndDown]
NonEEPath <- AllPaths[UpAndDown]

WithEEPath <- grep("EE",AllPaths)
DynPath <- AllPaths[-WithEEPath]

if(is.null(Paths)&OnlyDynamic==TRUE)PathsConsider <- DynPath
if(is.null(Paths)&OnlyDynamic==FALSE)PathsConsider <- NonEEPath
if(!is.null(Paths))PathsConsider <- Paths

HMMCall <- rownames(SigPPLarge)[which(SigPPLarge[,AllEEPath]<=FDR & SigMaxValLarge>cutoff
				  & SigMAPLargeChar%in%PathsConsider)]

Mat <- cbind(SigMAPLargeChar[HMMCall],round(SigMaxValLarge[HMMCall],4))
rownames(Mat) <- HMMCall
colnames(Mat) <- c("Most_Likely_Path","Max_PP")

MatOrder1 <- order(Mat[,2],decreasing=TRUE)
Mat1 <- Mat[MatOrder1,]
if(length(Mat1)==0)stop("No DE genes identified under defined FDR")
List <- sapply(PathsConsider,function(i){
	    tt <- which(Mat1[,1]==i)
	    if(length(tt)>0){
	    t2 <- matrix(Mat1[tt,],ncol=2)
	    rownames(t2) <- rownames(Mat1)[tt]
	    t2
				  }})

Length <- sapply(List,function(i){if(length(i)>0)tt <- nrow(i)
	      else tt <- 0
	      tt})

Names <- sapply(List,function(i){if(length(i)>0)tt <- rownames(i)
	      else tt <- NULL
	      tt})

Out <- list(Overall=Mat1,EachPath=List, NumEach=Length, EachPathNames=Names )
}

