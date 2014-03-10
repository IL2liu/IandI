library(gplots)
library(affy)
load("./data/annotation.RData/OrthoAnn.RData")

################################################################################
#               Now combine all of the data and statistics                     #
################################################################################

ids0 = cbind( rep(names(OrthoAnn$mm$hs),times=sapply(OrthoAnn$mm$hs,length)), 
             as.character(unlist(OrthoAnn$mm$hs)) )

##### data sets

# preprocessed data sets
load("./DiffProbeSets/data/RData/esets.RData")
eset.mm = esets$all
load("./InVitro/DiffProbeSets/data/RData/bppsl.eset.RData")
eset.hs = bppsl.eset
ids = ids0[ ids0[,1] %in% featureNames(eset.mm) | ids0[,2] %in% featureNames(eset.hs), ]
nmm = ncol(eset.mm)
nhs = ncol(eset.hs)

# Make the data matrix
comb.dat = as.data.frame(matrix(NA, nrow=nrow(ids), ncol=2+nmm+nhs))
comb.dat[,1:2] = ids
colnames(comb.dat) = c("Mouse","Human",colnames(eset.mm),colnames(eset.hs))
namat = comb.dat

# mouse
yn = comb.dat[,1] %in% featureNames(eset.mm)
ids.mm = comb.dat[ yn, 1]
comb.dat[yn,2+1:nmm] = exprs(eset.mm[ids.mm, ])

# human
yn = comb.dat[,2] %in% featureNames(eset.hs)
ids.hs = comb.dat[ yn, 2]
comb.dat[yn,2+nmm+1:nhs] = exprs(eset.hs[ids.hs, ])

##### fold change

# combine the fold changes
load("./DiffProbeSets/data/RData/cyberT.results.RData")
lfc.mm = cyberT.results$logFC
load("./InVitro/DiffProbeSets/data/RData/cybert.hs.results.RData")
lfc.hs = cybert.hs.results$fold
nmm = ncol(lfc.mm)
nhs = ncol(lfc.hs)

comb.lfc = namat[,1:(2+nmm+nhs)]
colnames(comb.lfc) = c("Mouse","Human",colnames(lfc.mm),colnames(lfc.hs))
namat2 = comb.lfc

# mouse
yn = comb.lfc[,1] %in% rownames(lfc.mm)
ids.mm = comb.lfc[ yn, 1]
comb.lfc[yn,2+1:nmm] = lfc.mm[ids.mm, ]

# human
yn = comb.lfc[,2] %in% rownames(lfc.hs)
ids.hs = comb.lfc[ yn, 2]
comb.lfc[yn,2+nmm+1:nhs] = lfc.hs[ids.hs, ]


##### t statistics
load("./DiffProbeSets/data/RData/cyberT.results.RData")
tstat.mm = cyberT.results$tstats
load("./InVitro/DiffProbeSets/data/RData/cybert.hs.results.RData")
tstat.hs = cybert.hs.results$bayesT

comb.tstat = namat2
yn = comb.tstat[,1] %in% rownames(tstat.mm)
ids.mm = comb.tstat[ yn, 1]
comb.tstat[yn,2+1:nmm] = tstat.mm[ids.mm, ]

yn = comb.tstat[,2] %in% rownames(tstat.hs)
ids.hs = comb.tstat[ yn, 2]
comb.tstat[yn,2+nmm+1:nhs] = tstat.hs[ids.hs, ]


##### p values
load("./DiffProbeSets/data/RData/cyberT.results.RData")
pval.mm = cyberT.results$pval
load("./InVitro/DiffProbeSets/data/RData/cybert.hs.results.RData")
pval.hs = cybert.hs.results$pVal

comb.pval = namat2
yn = comb.pval[,1] %in% rownames(pval.mm)
ids.mm = comb.pval[ yn, 1]
comb.pval[yn,2+1:nmm] = pval.mm[ids.mm, ]

yn = comb.pval[,2] %in% rownames(pval.hs)
ids.hs = comb.pval[ yn, 2]
comb.pval[yn,2+nmm+1:nhs] = pval.hs[ids.hs, ]
comb = list( pvals = comb.pval, dat = comb.dat, lfc = comb.lfc, tstat = comb.tstat )


# remove(list=setdiff(ls(),c(beginls,"comb")))
save(comb, file="./data/RData/MmAndHsCombined.RData")





  


