load("./InVitro/data/RData/abatch.RData")
library(limma)
source("./InVitro/Z-ExploreData/pcaLibrary.R")
library(affy)
library(affyPLM)
#load("./InVitro/data/annotation.RData/CustomPS.RData")

################################################################################
#               Create CDF environment from annotation list                    #
################################################################################

# download the cdf from brain array
# install.packages("/home/kevin/Downloads/hgu133plus2hsensgcdf_17.0.0.tar.gz")
library(hgu133plus2hsensgcdf)

hs = as.list(hgu133plus2hsensgcdf)
names(hs) = gsub("_at","",names(hs))
hgu133plus2ENSG = list2env(hs)

abatch@cdfName = "hgu133plus2ENSG"

################################################################################
#                              Preprocessing                                   #
################################################################################

# Set up the design matrix for every probe
a = interaction( pData(abatch)$Toxin, pData(abatch)$Time )
r = factor(pData(abatch)$Run)
dd = model.matrix( ~0+a+r )
colnames(dd) = gsub("a","",colnames(dd))

# With loess normalization
abatch.bgc = bg.correct.rma(abatch)
abatch.norm = normalize.AffyBatch.loess(abatch.bgc,type="pmonly",verbose=TRUE)
eset.pss = rma(abatch.norm, background=FALSE, normalize=FALSE)
eset.psn = normalize.ExpressionSet.loess(eset.pss)
fit = lmFit(eset.psn,dd)

# batch correction from loess normalization
esetfit = eset.psn
newset = exprs(esetfit)
for( i in 1:nrow(exprs(esetfit))) {
  resid = exprs(esetfit)[i,] - fit$design %*% fit$coefficients[i,]
  newset[i,] =  fit$design[,1:9] %*% fit$coefficients[i,1:9] + resid
}
bppsl.eset = esetfit
exprs(bppsl.eset) = newset

pca.plot(bppsl.eset)
save(bppsl.eset, file="./InVitro/DiffProbeSets/data/RData/bppsl.eset.RData")


