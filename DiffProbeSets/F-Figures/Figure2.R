library(ggplot2)
library(reshape)

################################################################################
#                           Venn diagrams                                      #
################################################################################

load("./DiffProbeSets/data/RData/cyberT.results.RData")
library(gplots)
distfun = function(x) as.dist( 1 - cor(t(x)) )
hclustfun = function(x) hclust(x,method="ward")

qvals = apply(cyberT.results$pval,2,function(x) p.adjust(x,method="BH"))
dif = as.data.frame(qvals[,-6]<0.01)

venn(dif[,1:3]) # 2 hours
venn(dif[,4:5]) # 6 hours
venn(dif[,6:7]) # 16 hours
venn(dif[,c(1,4,6)]) # TcdA
venn(dif[,c(2,5,7)]) # TcdB


################################################################################
#                           Large heatmap                                      #
################################################################################

load("./data/RData/cyberT.results.RData")
library(gplots)
distfun = function(x) as.dist( 1 - cor(t(x)) )
hclustfun = function(x) hclust(x,method="ward")

qvals = apply(cyberT.results$pval,2,function(x) p.adjust(x,method="BH"))
diff.probes = names( which( apply( qvals[,-6]<0.01, 1, any ) ) )
#diff.probes = names( which( apply( qvals[,]<0.01, 1, any ) ) )

FC = cyberT.results$logFC[diff.probes,]
colnames(FC) = c("A2","B2","AB2","A6","B6","AB6","A16","B16")
FC = FC[,c("B2","A2","AB2","B6","A6","AB6","B16","A16")]
exp.mat = FC[ , ]

heatmap.2( exp.mat,
           col=colorpanel(100,"blue","white","red"),
           breaks=seq(-5,5,length.out=101),
           symbreaks=TRUE,
           trace="none",
           Colv=NULL,
           dendrogram="none",
           distfun=distfun,
           hclustfun=hclustfun,
           density.info="none",
           labRow="")

################################################################################
#                      Correlation coefficients                                #
################################################################################

load("./DiffProbeSets/data/RData/cyberT.results.RData")
library(gplots)

logFC = cyberT.results$logFC
colnames(logFC) = c("A2","B2","AB2","A6","B6","AB6","A16","B16")
logFC = logFC[,c("B2","A2","AB2","B6","A6","AB6","B16","A16")]
cors = cor(logFC)

heatmap.2(cors,col=bluered,symbreak=TRUE,Rowv=NULL,
          Colv=NULL,dendrogram="none",trace="none",
          cellnote=round(cors,2),breaks=1000,
          density.info="none",asp=1)










