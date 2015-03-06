#' @title Likelihood function of the Beta-Negative Binomial HMM Model
#' @usage LikefunNBHMM(ParamPool, InputPool)
#' @param ParamPool The parameters that will be estimated in EM.
#' @param InputPool The control parameters that will not be estimated in EM
#' @author Ning Leng
#' @examples
#' data(GeneExampleData)
#' tmp <- GeneExampleData[1:10,]
#' In <- list(tmp,1,5,10,3,tmp,rep(1,15),as.factor(rep(1:5,each=3)), 10,cbind(rep(.5,10),rep(1,10),rep(2,10)))
#' Start <- c(1,1)
#' LikefunNBHMM(Start,In)
#' @details The likelihood function of the Beta-Negative Binomial HMM model used in EBSeqHMM.
#' EBSeqHMM uses optim() function to obtain the optimal estimates that minimizes the likelihood.
#' @return optimal estimates of the parameters of interest
#######################
# Optim fun
#######################
	LikefunNBHMM <- function(ParamPool, InputPool){

	a1 <- proc.time()[3]
	DataListIn.unlist <- InputPool[[1]]
	NumGene <- InputPool[[4]]
	AlphaMat <- InputPool[[2]]
	NumPoints <- InputPool[[3]]
	NumTranStage <- InputPool[[5]]
	RPairsListIn.unlist <- InputPool[[6]] 
	sizeFactors <- InputPool[[7]]

  Conditions <- InputPool[[8]] 
	NumOfEachGroupIn <- InputPool[[9]]

	FCParam <- InputPool[[10]]
	Vect <- ParamPool
	AlphaIn <- Vect[1]
	BetaIn <- Vect[2:(1+length(NumOfEachGroupIn))]




	LogBraw <- sapply(1:NumPoints,function(pp)
			 sapply(1:NumTranStage,function(ss){
	idx1 <- which(Conditions==levels(Conditions)[pp])
	idx2 <- which(Conditions==levels(Conditions)[pp+1])	

	t1 <-  matrix(DataListIn.unlist[,idx1],nrow=dim(DataListIn.unlist)[1])
	t2 <-  matrix(DataListIn.unlist[,idx2],nrow=dim(DataListIn.unlist)[1])
	tmpdata <- cbind(t1,t2)
	rownames(tmpdata) <- rownames(DataListIn.unlist)
	Rtmp1 <- outer(RPairsListIn.unlist[,pp],sizeFactors[idx1])			
	Rtmp2 <- outer(RPairsListIn.unlist[,pp],sizeFactors[idx2])			
	Log <- f0(tmpdata,
					        AlphaIn=AlphaIn, BetaIn=BetaIn, 
									cbind(Rtmp1*FCParam[ss,1],
										Rtmp2*FCParam[ss,2]), NumOfGroups=NumOfEachGroupIn, log=TRUE)
	
}),simplify=FALSE)	

	LogBraw2 <- do.call(rbind,LogBraw)
	LogB <- as.vector(LogBraw2)
	a2 <- proc.time()[3]
#cat(a2-a1)
#cat(" ")
	
	Sumgs <- LogB*AlphaMat
#	-sum(Sumgs)
#Use=which(Sumgs< quantile(Sumgs,.95,na.rm=TRUE)  & Sumgs> quantile(Sumgs,0,na.rm=TRUE), arr.ind=TRUE)
Use <- which(!is.na(Sumgs) & abs(Sumgs)!=Inf, arr.ind=TRUE)
	-mean(Sumgs[Use])
}

