#' @title Obtain DE gene/isoform list at a certain FDR 
#' @usage GetDECalls(EBSeqHMMOut,FDR=.05)
#' @param EBSeqHMMOut output from EBSeqHMMTest function
#' @param FDR Target FDR; default is 0.05 
#' @note Output: output a list of genes that are DE in at least one condition in an RNA-seq
#' experiment with multiple ordered conditions
#' @author Ning Leng
#' @examples 
#' data(GeneExampleData)
#' CondVector <- rep(paste("t",1:5,sep=""),each=3)
#' Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
#' Sizes <- MedianNorm(GeneExampleData)
#' EBSeqHMMGeneOut <- EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
#'           UpdateRd=2)
#' GeneDECalls <- GetDECalls(EBSeqHMMGeneOut, FDR=.05)
#' @details Function GetDECalls() can be used to obtain a list of DE genes/isoforms
#' with user specific cutoffs. To obtain a list of DE genes/isoforms with
#' a target FDR alpha, the user may specify FDR=alpha. 
#' @return a list of genes/isoforms that are identified as DE under the target FDR, shown are their names and PPs;




GetDECalls <- function(EBSeqHMMOut, FDR=.05){

SigPPLarge <- EBSeqHMMOut$MgAllPP # PP matrix
SigMAPLargeChar <- EBSeqHMMOut$MgAllMAPChar # 
SigMaxValLarge <- EBSeqHMMOut$MgAllMaxVal

AllPaths <- colnames(EBSeqHMMOut$MgAllPP)
WithUp <- grep("Up",AllPaths)
WithDown <- grep("Down",AllPaths)
UpAndDown <- union(WithUp, WithDown)
AllEEPath <- AllPaths[-UpAndDown]
HMMNonEECall <- rownames(SigPPLarge)[which(SigPPLarge[,AllEEPath]<=FDR)]

Mat <- cbind(SigMAPLargeChar[HMMNonEECall],round(SigMaxValLarge[HMMNonEECall],4))
rownames(Mat) <- HMMNonEECall
colnames(Mat) <- c("Most_Likely_Path","Max_PP")

MatOrder1 <- order(Mat[,2],decreasing=TRUE)
Mat1 <- Mat[MatOrder1,]

Out <- Mat1

}

