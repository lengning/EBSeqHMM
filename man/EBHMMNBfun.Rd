\name{EBHMMNBfun}
\alias{EBHMMNBfun}
\title{Baum-Welch algorithm for a single hidden markov chain}
\usage{
EBHMMNBfun(Data,NgVector=NULL,Conditions, sizeFactors,
	PriorFC=1.5,homo=TRUE, maxround=5,
	Pi0=NULL, Tran=NULL,NoTrend=FALSE, NumTranStage=3,
	FCParam=NULL, AlphaIn=NULL,BetaIn=NULL,
	StateNames=c("Up","NC","Down"),
	EM=TRUE, UpdateParam=TRUE, Print=TRUE,
	OnlyQ=FALSE,WithinCondR=TRUE,
	PenalizeLowMed=TRUE, PenalizeLowMedQt=.2,PenalizeLowMedVal=10)
}
\arguments{
  \item{Data}{input data, rows are genes/isoforms and
  columns are samples}

  \item{NgVector}{Ng vector; NULL for gene level data}

  \item{Conditions}{A factor indicates the condition
  (time/spatial point) which each sample belongs to.}

  \item{sizeFactors}{a vector indicates library size
  factors}

  \item{Tran}{initial values for transition matrices}

  \item{Pi0}{initial values for starting probabilities}

  \item{NumTranStage}{number of states}

  \item{PriorFC}{target FC for gridient change}

  \item{StateNames}{name of the hidden states}

  \item{homo}{whether the chain is assumed to be
  homogenious}

  \item{maxround}{max number of iteration}

  \item{AlphaIn,BetaIn}{If the parameters are known and the
  user doesn't want to estimate them from the data, user
  may specify them here.}

  \item{NoTrend}{if NoTrend=TRUE, initial transition
  probabilities will be set to be equal}

  \item{FCParam}{not in use}

  \item{EM}{Whether estimate the prior parameters of the
  beta distribution by EM}

  \item{UpdateParam}{Whether update starting probabilities
  and transition probabilities}

  \item{OnlyQ}{If OnlyQ=TRUE, the function will only return
  estimated q's.}

  \item{WithinCondR}{By defining WithinCondR=TRUE,
  estimation of r's are obtained within each condition.
  (WithinCondR=FALSE is not suggested here)}

  \item{Print}{Whether print the elapsed-time while running
  the test.}

  \item{PenalizeLowMed,PenalizeLowMedQt,PenalizeLowMedVal}{Transcripts
  with median quantile < = PenalizeLowMedQt will be
  penalized}
}
\value{
MAPTerm: the most likely path of each gene/isoform.
MAPTermNum: numeric version of MAPTerm.

AllTerm: all possible expression paths considered in the
model. PP: posterior probability of being each expression
path.

WhichMax: index of the most likely path. Allf: prior
probability of being each path.

Pi0Track: estimated starting probabilities of each
iteration.

TranTrack: estimated transition probabilities of each
iteration.

AlphaTrack, BetaTrack: estimated alpha and beta(s).

LLAll=PostSumForLL.Sum: log likelihood of the model.
}
\description{
Baum-Welch algorithm for a single hidden markov chain
}
\details{
EBHMMNBfun() function implements the Balm-Welch algorithm
that estimates the starting probabilities and transition
probabilities of a single hidden Markov model. Here the
emission distribution of each gene is assumed to be a
Beta-Negative Binomial distribution with parameters (r_g,
alpha, beta) , in which alpha and beta are shared by all
the genes and r_g is gene specific. If not specified, r_g,
alpha and beta will be estimated using method of moments.
For isoform data, we assume isoforms from the same Ig group
share the same beta^Ig. alpha is shared by all the isoforms
and r_gi is isoform specific. The user also needs to
specify an expected FC.
}
\examples{
data(GeneExampleData)
CondVector <- rep(paste("t",1:5,sep=""),each=3)
Conditions <- factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
Sizes <- MedianNorm(GeneExampleData)
tmp <- EBHMMNBfun(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
          maxround=2, OnlyQ=TRUE)
}
\author{
Ning Leng
}

