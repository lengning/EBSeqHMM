#' @title Calculate the prior predictive distribution of the Beta-Negative Binomial  model
#' @usage f0(Input, AlphaIn, BetaIn, EmpiricalR, NumOfGroups, log)
#' @param Input expression values
#' @param AlphaIn,BetaIn,EmpiricalR The parameters estimated from last
#'          iteration of EM.
#' @param NumOfGroups How many transcripts within each Ng group
#' @param log If set as TRUE, the output will in log scale.
#' @author Ning Leng
#' @examples
#'     f0(matrix(rnorm(100,100,1),ncol=10), .5, .6,
#'            matrix(rnorm(100,200,1),ncol=10), 100, TRUE)
#' @details Function f0() will calculate the Beta-Negative Binomial prior predictive probability for a given set of parameters
#' @return output a numeric vector, each element shows the prior predictive probability of one gene/isoform

f0 <-
function(Input, AlphaIn, BetaIn, EmpiricalR, NumOfGroups, log)
{	
		BetaVect <- do.call(c,sapply(1:length(BetaIn),function(i)rep(BetaIn[i],NumOfGroups[i]),simplify=FALSE))
		SampleNum <- dim(Input)[2]
		#Product part
		ChooseParam1 <- round(Input+EmpiricalR-1)
		roundInput <- round(Input)
		EachChoose0 <- matrix(sapply(1:SampleNum, function(i)lchoose(ChooseParam1[,i], roundInput[,i])),ncol=SampleNum)
	
	# numerical approximation to rescue -Inf ones
    NoNegInfMin <- min(EachChoose0[which(EachChoose0!=-Inf)])
    NoPosInfMax <- max(EachChoose0[which(EachChoose0!=Inf)])
    EachChoose <- EachChoose0
    EachChoose[which(EachChoose0==-Inf, arr.ind=TRUE)] <- NoNegInfMin
		EachChoose[which(EachChoose0==Inf, arr.ind=TRUE)] <- NoPosInfMax
	
		
		SumEachIso <- rowSums(Input)
		param1 <- AlphaIn + rowSums(EmpiricalR)
		param2 <- BetaVect + SumEachIso
		LogConst <- rowSums(EachChoose)+lbeta(param1, param2)-lbeta(AlphaIn, BetaVect)


		if (log==FALSE) FinalResult <- exp(LogConst)
		if (log==TRUE) FinalResult <- LogConst
    FinalResult
}

