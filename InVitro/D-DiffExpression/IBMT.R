# load the batch processed data
library(limma)
source("./InVitro/D-DiffExpression/ibmtFunction.R")
source("./C-Annotations/CollapseIDs.R")
load("./InVitro/data/RData/bppsl.eset.RData")
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
save(ibmt.results,file="./InVitro/data/RData/ibmtResults.RData")
ibmt.toptable(ibmt.results, 1)

################################################################################
#                            Collapsed IBMT                                    #
################################################################################

# collapse the complete data set
collapsed.ibmt.results = collapsedIBMT( eset, HsAnn$Symbol$Affy, dd, cm )
save(collapsed.ibmt.results,file="./InVitro/data/RData/collapsed.ibmtResults.RData")

res = collapsed.ibmt.results
ibmt.toptable(res, 2)


################################################################################
#                       multiContrast Collapsed IBMT                           #
################################################################################

# Separate the data by contrasts
mc.ibmt = vector("list",ncol(cm))

# run IBMT for each collapsed contrast
for( i in 1:ncol(cm) ) {
  arrays = apply( dd[ , cm[,i]!=0 ]!=0, 1, any )
  subeset = eset[,arrays]
  subdd = dd[arrays,cm[,i]!=0]
  subcm = t(t(cm[cm[,i]!=0,i]))
  mc.ibmt[[i]] = collapsedIBMT( subeset, HsAnn$Symbol$Affy, subdd, subcm )
}

# now combine the multiple lists as much as possible
combined = list()

# combine all the coefficients
aa = sapply(mc.ibmt, function(x) x$coefficients)
rownames(aa) = rownames(mc.ibmt[[1]]$coefficients)
colnames(aa) = colnames(cm)
combined$coefficients = aa

# combine all the other vectors
classes = sapply(mc.ibmt[[1]], class)
for( i in which(classes=="integer" | classes=="numeric") ) {
  out = matrix(sapply( mc.ibmt, function(x) x[[i]] ),ncol=ncol(cm))
  rownames(out) = names(mc.ibmt[[1]][[i]])
  colnames(out) = colnames(cm)
  combined[[  names(mc.ibmt[[1]])[i]  ]]  = out
}
mc.ibmt$combined = combined

save(mc.ibmt,file="./InVitro/data/RData/mc.ibmt.RData")


