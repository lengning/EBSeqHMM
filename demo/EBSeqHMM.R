##################
# Gene level analysis
##################
data(GeneExampleData)
str(GeneExampleData)

CondVector=rep(paste("t",1:5,sep=""),each=3)
print(CondVector)
Conditions=factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
str(Conditions)


Sizes=MedianNorm(GeneExampleData)

EBSeqHMMGeneOut=EBSeqHMMTest(Data=GeneExampleData, sizeFactors=Sizes, Conditions=Conditions,
														            UpdateRd=5)

# Identify DE genes
GeneDECalls=GetDECalls(EBSeqHMMGeneOut, FDR=.05)
head(GeneDECalls)
str(GeneDECalls)

# Classify DE genes into expression paths
GeneConfCalls=GetConfidentCalls(EBSeqHMMGeneOut, FDR=.05,cutoff=.5, OnlyDynamic=T)
print(GeneConfCalls$EachPath[1:4])

Path4=GeneConfCalls$EachPath[["Down-Down-Up-Up"]]
print(Path4)

head(GeneConfCalls$Overall)
print(GeneConfCalls$NumEach)
str(GeneConfCalls$EachPathGeneNames)


##################
# Isoform level analysis
##################
data(IsoExampleList)
str(IsoExampleList)
IsoExampleData=IsoExampleList$IsoExampleData
str(IsoExampleData)
IsoNames=IsoExampleList$IsoNames
IsosGeneNames=IsoExampleList$IsosGeneNames

IsoSizes=MedianNorm(IsoExampleData)

NgList=GetNg(IsoNames, IsosGeneNames)
IsoNgTrun=NgList$IsoformNgTrun
IsoNgTrun[c(1:3,101:103,161:163)]

CondVector=rep(paste("t",1:5,sep=""),each=3)
Conditions=factor(CondVector, levels=c("t1","t2","t3","t4","t5"))
str(Conditions)

EBSeqHMMIsoOut=EBSeqHMMTest(Data=IsoExampleData,
									          NgVector=IsoNgTrun,
									          sizeFactors=IsoSizes, Conditions=Conditions,
									          UpdateRd=5)

# Identify DE isoforms
IsoDECalls=GetDECalls(EBSeqHMMIsoOut, FDR=.05)
str(IsoDECalls)
head(IsoDECalls)

# Classify DE isoforms into expression paths
IsoConfCalls=GetConfidentCalls(EBSeqHMMIsoOut, FDR=.05)
head(IsoConfCalls$Overall)
str(IsoConfCalls$EachPath[1:4])
str(IsoConfCalls$EachPathGeneNames)
