# Compare the proteomics study to our data

aa = read.csv("./data/pr300973q_si_001.csv",skip=2)
bb = aa[, c(1,4,6,7:15,31,32,34)]

# average all of the log2 scores
cc = bb[,4:6]
dd = apply(cc,2,as.numeric) # rTcdA over control

# calculate fold changes
foo = function(x) {
  num.nas = sum(is.na(x))
  num.zeros = sum(x==0,na.rm=TRUE)
  if( num.nas==2 & num.zeros==1 || num.nas==3 ) { # if 2 NAs and zero OR all NAs, make the average NA
    return( NA )
  } else {
    return( mean(x,na.rm=TRUE) ) # remove the NA values and average
  }
}
fc.av = apply(dd,1,foo)
#fc.av = rowMeans(dd,na.rm=TRUE)

# log folc change
lfc = log2(fc.av)
syms = as.character(bb$Gene.Names)

# clean up the data
rr = data.frame(syms,lfc,stringsAsFactors=FALSE)
rr2 = rr[ rr$syms != "", ]

################################################################################
#                Compare Caco2 proteome to HCT8 transcription                  #
################################################################################

# now worry about the mappings from protein to gene
load("./InVitro/DiffProbeSets/data/RData/cybert.hs.collapsed.results.RData")
genelfc = cybert.hs.collapsed.results$fold
gene.syms = rownames(genelfc)
symsplit = strsplit(rr2$syms,";")
foo1 = function(x) {
  indata = x %in% gene.syms
  if( all(!indata) ) {
    return("Nothing")
  } else {
    return(x[indata])
  }
}
symsplit2 = lapply(symsplit,foo1)
table(sapply(symsplit2,length))

# first just look at correlation with genes that map 1 to one
symsplit3 = sapply(symsplit2,function(x) paste(x,collapse=";"))
rr3 = data.frame( symsplit3, rr2$lfc )
rr4 = rr3[ rr3$symsplit3!="Nothing", ]

# combine the data points
nn = cbind( genelfc[ rr4$symsplit3, ], rr4$rr2.lfc )
oo = na.omit(nn)
cor(oo, method="spearman")
oo.hum = oo
  
# visualize the correlations
#plot(oo[,3],oo[,7])
library(gplots)
heatmap.2(cor(oo),col=colorpanel(100,"blue","white","red"),symbreaks=TRUE,
          trace="none",cellnote=round(cor(oo,method="spearman"),2),
          Colv=FALSE)


################################################################################
#               Compare Caco2 proteome to in vivo transcription                #
################################################################################

load("./DiffProbeSets/data/RData/cyberT.results.RData")
load("./data/annotation.RData/MsAnn.RData")
load("./data/annotation.RData/OrthoAnn.RData")
load("./InVitro/data/annotation.RData/HumanAnnotations.RData")
source("./C-Annotations/TransitiveMapping.R")
source("./C-Annotations/CollapseIDs.R")
maap = TransitiveMapping( HsAnn$Ensembl$Symbol, OrthoAnn$hs$mm  )
foo2 = function(x) sum(abs(x))
genelfc = collapseExprMatrix(cyberT.results$logFC, maap, foo2)
gene.syms = rownames(genelfc)

symsplit = strsplit(rr2$syms,";")
foo1 = function(x) {
  indata = x %in% gene.syms
  if( all(!indata) ) {
    return("Nothing")
  } else {
    return(x[indata])
  }
}
symsplit2 = lapply(symsplit,foo1)
symsplit3 = sapply(symsplit2,function(x) paste(x,collapse=";"))
rr3 = data.frame( symsplit3, rr2$lfc )
rr4 = rr3[ rr3$symsplit3!="Nothing", ]

nn = cbind( genelfc[ rr4$symsplit3, ], rr4$rr2.lfc )
oo = na.omit(nn)
cor(oo, method="spearman")
oo.m = oo

plot(oo[,8],oo[,9])
library(gplots)
heatmap.2(cor(oo),col=colorpanel(100,"blue","white","red"),symbreaks=TRUE,
          trace="none",cellnote=round(cor(oo,method="spearman"),2))

################################################################################
#                   Combine all three data sets                                #
################################################################################

aa1 = oo.m[ rownames(oo.m) %in% rownames(oo.hum), ]
aa2 = cbind( aa1[,-ncol(aa1)], oo.hum[ rownames(aa1), ] )
heatmap.2(cor(aa2,method="spearman"),col=colorpanel(100,"blue","white","red"),
          symbreaks=TRUE,Colv=FALSE,Rowv=FALSE,dendrogram="none",
          trace="none",cellnote=round(cor(aa2,method="spearman")^2,2),
          asp=1)



