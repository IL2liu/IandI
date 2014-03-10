# A script that runs cyberT for all the different two class comparisons
load("./DiffProbeSets/data/RData/esets.RData")
load("./data/annotation.RData/MsAnnGO.RData")
load("./data/annotation.RData/MsAnn.RData")
load("./data/annotation.RData/GSAnn.RData")
source("./E-Enrichment/cactus.R")

# Arrays were preprocessed in groups designated by time point (e.g. 2h, 6h or 16h)
# Here combine them only for the sake of simpler code
eset = new("ExpressionSet")
e1 = esets[["2hr"]]
e2 = esets[["6hr"]]
e3 = esets[["16hr"]]
exprs(eset) = cbind(exprs(e1),exprs(e2),exprs(e3))
pData(eset) = rbind(pData(e1),pData(e2),pData(e3))

comparisons = list( c("A2","Sham2"), c("B2","Sham2"), c("AB2","Sham2"),
                    c("A6","Sham6"), c("B6","Sham6"),
                    c("A16","Sham16"), c("B16","Sham16") )
mapping=MsAnn$MGI$Ensembl

# Run camera on different groups of gene sets
result = list()
result$mf = cactus(eset, MsAnnGO$MF$MGI, comparisons, mapping )
result$bp = cactus(eset, MsAnnGO$BP$MGI, comparisons, mapping )
result$cc = cactus(eset, MsAnnGO$CC$MGI, comparisons, mapping )
result$go = cactus(eset, c(MsAnnGO$CC$MGI,MsAnnGO$BP$MGI,MsAnnGO$MF$MGI), 
                             comparisons, mapping )

cactus.result = result
save(cactus.result, file="./DiffProbeSets/data/RData/cactus.RData")

#i=1
#aa = result$go[,i]
#hist(aa)
#Term(names(sort(aa)[1:27]))

