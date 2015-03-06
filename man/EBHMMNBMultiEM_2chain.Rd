\name{EBHMMNBMultiEM_2chain}
\alias{EBHMMNBMultiEM_2chain}
\title{Run EBSeqHMM model with a fixed expected FC}
\usage{
EBHMMNBMultiEM_2chain(Data,
	NgVector=NULL, Conditions, AllTran=NULL,
	AllPi0=NULL, Terms=NULL,
	sizeFactors, NumTranStage=c(3,2),PriorFC=2,
	StateNames=c("Up","Down"),homo=FALSE,
	UpdateRd=5, PIBound=.9, UpdatePI=FALSE,Print=FALSE,
	WithinCondR=TRUE,
	PenalizeLowMed=TRUE, PenalizeLowMedQt=.1,PenalizeLowMedVal=10)
}
\arguments{
  \item{Data}{input data, rows are genes and columns are
  samples}

  \item{NgVector}{Ng vector; NULL for gene level data}

  \item{Conditions}{A factor indicates the condition
  (time/spatial point) which each sample belongs to.}

  \item{AllTran}{initial values for transition matrices}

  \item{AllPi0}{initial values for starting probabilities}

  \item{Terms}{Terms}

  \item{sizeFactors}{a vector indicates library size
  factors}

  \item{StateNames}{names of the hidden states}

  \item{NumTranStage}{number of states in two chains}

  \item{PriorFC}{target FC for gridient change}

  \item{homo}{whether the chain is assumed to be
  homogenious}

  \item{UpdateRd}{max number of iteration}

  \item{UpdatePI}{whether update the mixture proportion of
  two chains}

  \item{PIBound}{upper bound of the mixture proportion of
  the two chains}

  \item{Print}{Whether print the elapsed-time while running
  the test.}

  \item{WithinCondR}{By defining WithinCondR=TRUE,
  estimation of r's are obtained within each condition.
  (WithinCondR=FALSE is not suggested here)}

  \item{PenalizeLowMed,PenalizeLowMedQt,PenalizeLowMedVal}{Transcripts
  with median quantile < = PenalizeLowMedQt will be
  penalized}
}
\value{
Pi0Out: estimated starting probabilities of each iteration.

TranOut: estimated transition probabilities of each
iteration.

Pi: estimated probability of being each chain.

Alpha, Beta: estimated alpha and beta(s).

LLSum: log likelihood of the model.

QList: estimated q's.

MgAllPP: marginal PP for all paths.

MgAllMAPChar: Most likely path based on MgAllPP.

MgAllMaxVal: Highest PP based on MgAllPP.

PPMatW: Posterior probabilities of being each of the
chains.
}
\description{
Run EBSeqHMM model with a fixed expected FC
}
\details{
EBHMMNBMultiEM_2chain() function implements the EBSeqHMM
model to perform statistical analysis in an RNA-seq
experiment with ordered conditions. EBHMMNBMultiEM_2chain()
calls EBHMMNBfunForMulti() function to perform Balm-Welch
algorithm that estimates the starting probabilities and
transition probabilities. Here the emission distribution of
each gene is assumed to be a Beta-Negative Binomial
distribution with parameters (r_g, alpha, beta) , in which
alpha and beta are shared by all the genes and r_g is gene
specific. If not specified, r_g, alpha and beta will be
estimated using method of moments. For isoform data, we
assume isoforms from the same Ig group share the same
beta^Ig. alpha is shared by all the isoforms and r_gi is
isoform specific. The user also needs to specify an
expected FC. Function EBSeqHMMTest() runs several models
with varying FCs and returns the model with maximum
likelihood.
}
\examples{
data(GeneExampleData)
CondVector <- rep(paste("t",1:5,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
Sizes <- MedianNorm(GeneExampleData)
tmp <- EBHMMNBMultiEM_2chain(Data=GeneExampleData,sizeFactors=Sizes, Conditions=Conditions,
          UpdateRd=2)
}
\author{
Ning Leng
}

