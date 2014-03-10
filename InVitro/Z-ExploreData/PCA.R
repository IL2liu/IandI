load("./InVitro/data/RData/esets.RData")
#load("./InVitro/data/RData/esets.rma.RData")
library(affy)
library(affyPLM)
source("./InVitro/Z-ExploreData/pcaLibrary.R")

################ all samples

pca.plot( esets[[1]], 1, 2 )
pca.plot( esets[[1]], 2, 4 )

pca.plot( esets[[2]], 1, 2 )
pca.plot( esets[[2]], 2, 4 )

pca.plot( esets[[3]], 1, 2 )

pca.plot( esets$all, 1, 2 )
pca.plot( esets$all, 1, 3 )

pData(esets$all)[1:2,3] = 1:2


############# fold change
source("parseCyberT.R")
all.eset = esets.processed$all
ab6.eset = all.eset[ , pData( all.eset )$Groups == "AB6" ]
c6.eset = all.eset[ , pData( all.eset )$Groups == "C6" ]
ab6.fold = exprs(ab6.eset) - rowMeans(exprs(c6.eset))
all.fold = cbind(fold[,1:5],ab6.fold,fold[,6:7])
colnames(all.fold)[6] = "AB6-C6"

devSVG(file="../Figures/Figure1/PCA.fold.svg")
pca.fold( all.fold, 1, 2 )
dev.off()
