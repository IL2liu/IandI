load("./InVitro/data/RData/abatch.RData")
library(limma)
source("./InVitro/Z-ExploreData/pcaLibrary.R")
library(affy)
library(affyPLM)

# Set up the design matrix for every probe
a = interaction( pData(abatch)$Toxin, pData(abatch)$Time )
r = factor(pData(abatch)$Run)
dd = model.matrix( ~0+a+r )


################################################################################
#            Perform batch correction after probe set summarization            #
################################################################################

# With RMA
eset.rma = rma(abatch)
pca.plot(eset.rma)

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
save(bppsl.eset, file="./InVitro/data/RData/bppsl.eset.RData")


################################################################################
#              Before and after batch correction - PCA plots                   #
################################################################################

eset.rma = rma(abatch)
load("./InVitro/data/RData/bppsl.eset.RData")
eset.bc = bppsl.eset

# get pca results from each eset
pcaa = function( eset ) {
  centered = exprs(eset) - rowMeans(exprs(eset))
  selected = rowSums( centered == 0 ) != ncol(eset)
  eset.pca = eset[selected,]
  pca = prcomp( t(exprs(eset.pca)), center=TRUE, scale=TRUE )
  pca.info = summary(pca)
  pca.info
}
pca1 = pcaa(eset.rma)
pca2 = pcaa(eset.bc)
x1 = pca1$x[,1]
y1 = pca1$x[,2]
x2 = pca2$x[,1]
y2 = pca2$x[,2]
mraw = cbind( pData(eset.rma), x1, y1, x2, y2 )
mraw$Run = factor(mraw$Run)
mraw$Time = factor(mraw$Time)
mraw$group = gsub("[.]","",interaction(mraw$Toxin,mraw$Time))


library(ggplot2)
library(gridExtra)
p1 = ggplot(data=mraw) +
  geom_text(aes(label=group,x=x1,y=y1,color=Run))
p2 = ggplot(data=mraw) +
  geom_text(aes(label=group,x=x2,y=y2,color=Run))
grid.arrange(p1,p2,nrow=1)







