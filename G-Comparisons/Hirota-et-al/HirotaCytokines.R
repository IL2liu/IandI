############# Hirota et al paper comparison
library(reshape)
library(ggplot2)
library(gplots)
source("./G-Comparisons/Hirota-et-al/library.R")

################################################################################
#                   Parsing and formatting the data                            #
################################################################################

serum = getSerumData()
hirota = getHirotaData()
cells = getCellData()

################################################################################
#               Hirota cytokines versus cell & serum cytokines                 #
################################################################################

prots = intersect( colnames(serum$lfc), colnames(hirota$lfc))
allprots = t(rbind( serum$lfc[,prots], cells$lfc[,prots], hirota=hirota$lfc[,prots] ))

pairs(~.,allprots)
corplot(allprots)

################################################################################
#                 Hirota cytokines versus gene expression                      #
################################################################################

prots = rownames(na.omit(t(hirota$lfc)))
genes = getGeneData(prots)
hirota.v.genes = cbind( hirota=hirota$lfc[,prots], genes )

pairs(~.,hirota.v.genes)
corplot(hirota.v.genes)


################################################################################
#              Cell and serum cytokines versus gene expression                 #
################################################################################

prots = colnames(serum$lfc)
genes = getGeneData(prots)
cs.genes = cbind( t(rbind( serum$lfc[,prots], cells$lfc[,prots] )),genes )

pairs(~.,cs.genes)
corplot(cs.genes)

################################################################################
#                    Plotting the key comparisons                              #
#################################################################################

# For cells, Ccl4 not plotted. near detection level. high variance
cells$data

# Some genes that have expression levels very near detection levels.
# Explore the data
load("./DiffProbeSets/data/RData/esets.RData")
mgis = unlist( MsAnn$Symbol$MGI[ rownames(hirota.v.genes) ] )
arrays = pData(esets$all)$Toxin %in% c("A","Sham","AB") & 
          pData(esets$all)$Time ==6
raw = exprs(esets$all[mgis,arrays])
colnames(raw) = apply( pData(esets$all)[arrays,], 1, function(x) paste(x,collapse="") )
data = cbind( raw, cyberT.results$logFC[mgis,c(4,6)], cyberT.results$pval[mgis,c(4,6)] )
rownames(data) = names(mgis)
data

# custom plotting function
cplot = function(i,j,data=fig.cks,ylim=c(0,9),xlim=c(-2,4)) {
  plot(data[,i],data[,j],
       xlab=colnames(data)[i],
       ylab=colnames(data)[j],
       xlim=xlim, ylim=ylim, col="white",#bty="n",
       main=cor(data[,i],data[,j]))
  text(data[,i],data[,j],rownames(data))
  grid()
  abline(0,1)
  abline(h=0)
  abline(v=0)
}

# plot all of the comparisons
par(mfrow=c(2,2))
cplot(1,2,cs.genes[-7,c(3,8)],xlim=c(-0.1,3.2),ylim=c(-1,4.2)) 
cplot(2,1,hirota.v.genes[,c(1,7)],xlim=c(-0.3,4),ylim=c(-0.3,9.5))
cplot(1,2,allprots[-6,c(3,5)],xlim=c(-0.1,3.2),ylim=c(0,8.4))
cplot(1,2,allprots[,c(2,5)],xlim=c(-0.2,6.5),ylim=c(-0.2,8.4))






