# load the batch processed data
library(limma)
source("./InVitro/D-DiffExpression/ibmtFunction.R")
source("./C-Annotations/CollapseIDs.R")
load("./InVitro/DiffProbeSets/data/RData/bppsl.eset.RData")
load("./InVitro/data/annotation.RData/HumanAnnotations.RData")
eset = bppsl.eset

# Set up the design matrix
a = interaction( pData(eset)$Toxin, pData(eset)$Time )
dd = model.matrix( ~0+a )
colnames(dd) = gsub("a","",colnames(dd))
colnames(dd) = gsub("[.]","",colnames(dd))
cm = makeContrasts( A2-C2, B2-C2, A6-C6, B6-C6, A24-C24, B24-C24, levels=dd )

################################################################################
#                             Regular IBMT                                     #
################################################################################

# fit linear models to the data
fit = lmFit(eset, dd)
fit2 = contrasts.fit(fit,cm)

# perform IBMT for differential expression
ibmt.results = IBMT( fit2, 1:ncol(fit2) )
save(ibmt.results,file="./InVitro/DiffProbeSets/data/RData/ibmtResults.RData")
ibmt.toptable(ibmt.results, 1, probeann="ENSG")

################################################################################
#                   Collapsed to Gene Symbol, IBMT                             #
################################################################################

# collapse the complete data set
collapsed.ibmt.results = collapsedIBMT( eset, HsAnn$Symbol$Ensembl, dd, cm )
save(collapsed.ibmt.results,file="./InVitro/DiffProbeSets/data/RData/collapsed.ibmtResults.RData")