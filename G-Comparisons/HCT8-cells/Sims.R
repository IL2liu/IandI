# Identify genes which are similar and then discuss which ones
# they are

load("./data/RData/MmAndHsCombined.RData")
load("./data/annotation.RData/MsAnn.RData")
source("./C-Annotations/TransitiveMapping.R")
map = TransitiveMapping( MsAnn$MGI$Ensembl, MsAnn$MGI$Symbol )
library(gplots)

# which rows belong to which data set
mm = 3:10
hs = 11:16
mmsyms = comb$lfc[,1]
hssyms = comb$lfc[,2]

# separate mouse and human data
mmd = comb$lfc[,mm]
hsd = comb$lfc[,hs]
alld = comb$lfc[,c(mm,hs)]

################################################################################
#         Choose by high differential expression in both data sets             #
################################################################################

# find which genes are overall commonly regulated across all time points
mmm = rowMeans(abs(mmd))
hsm = rowMeans(abs(hsd))

# which genes are regulated in both data sets
reg = mmm>4.4 & hsm>4.4
reg = mmm>log2(1.5) & hsm>log2(1.5)
same = sign(rowMeans(mmd)) == sign(rowMeans(hsd))
genes = which( reg & same )

# what are the gene symbbols
ensgs = mmsyms[genes]
syms = sapply( map[ensgs], function(x) paste(x,collapse=";") )

# plot the expression data
hm = data.matrix(alld[genes,])
rownames(hm) = syms
hm = hm[-13,] # very bad hack because Alpi occurs twice
heatmap.2(hm,col=colorpanel(100,"blue","white","red"),symbreaks=TRUE,
          trace="none", Colv=FALSE, dendrogram="none",
          cellnote=round(sign(hm)*2^abs(hm),1))

################################################################################
#                   Last, do a correlation heatmap                             #
################################################################################

# function to do correlation for data sets with NAs
na.cor = function( x ) {
  out = matrix(0,ncol(x),ncol(x))
  for( i in 1:ncol(x) ) {
    for( j in 1:ncol(x) ) {
      out[i,j] = cor(na.omit(cbind(x[,i],x[,j])),method="spearman")[2]
    }
  }
  colnames(out) = rownames(out) = colnames(x)
  return(out)
}

x = data.matrix(comb$lfc[,3:ncol(comb$lfc)])
cc = na.cor(x)
heatmap.2(cc,col=colorpanel(100,"blue","white","red"),symbreaks=TRUE,
          trace="none", Colv=FALSE, Rowv=FALSE, dendrogram="none",
          cellnote=round(cc^2,2),key=FALSE,asp=1)








