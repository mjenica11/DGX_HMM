# Load libraries
library(EBSeq)
library(DESeq2)

# Load simulated data
data(GeneMat)

# Data is matrix of gene counts; genes are rows, cols are conditions
str(GeneMat)

# EBSeq requires library size factors to adjust for sequencing depth among samples
Sizes=MedianNorm(GeneMat)

# Simulate five samples to be in condition 1 and condition 2
# Detect DE genes
# Number of iterations of EM algorithm set by maxround
EBOut <- EBTest(Data=GeneMat, Conditions=as.factor(rep(c("C1", "C2"), each=5)), sizeFactors=Sizes, maxround=5)

# Get list of DE genes and posterior probs of being DE; FDR = 'Fasle Discovery Rate'
EBDERes <- GetDEResults(EBOut, FDR=0.05)
str(EBDERes$DEfound)

# PPEE/PPDE col; posterior prob of being EE or DE 
head(EBDERes$PPMat)
str(EBDERes$Status)

# Calculate FC (fold change); plot posterior vs raw FC
GeneFC=PostFC(EBOut)
str(GeneFC)
PlotPostVsRawFC(EBOut,GeneFC)

# Check convergence; assumed prior dist of qg^c is Beta(alpha,beta)
# EM estimates the hyper-parameters alpha,beta and the mixture parameter p
EBOut$Alpha
EBOut$Beta
EBOut$P

# Check model fit; empirical q vs simulated q's from Beta prior distribution
par(mfrow=c(1,2))
QQP(EBOut)

# Density plot of " " " from Beta prior dist
par(mfrow=c(1,2))
DenNHist(EBOut)






