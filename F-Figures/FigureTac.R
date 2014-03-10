library(gplots)
library(ggplot2)
library(reshape)
library(gridExtra)
.pt = 0.47 # found this by trial and error to work on my machine
.ptText = 3/8.5
source("./C-Annotations/CollapseIDs.R")
source("./C-Annotations/TransitiveMapping.R")
load("./data/annotation.RData/MsAnn.RData")
load("./data/annotation.RData/MsAnnGO.RData")
load("./data/annotation.RData/GSAnn.RData")
load("./data/RData/cyberT.results.RData")
library(GO.db)

genes = c("Tac1","Tac2","Tac4","Tacr1","Tacr2","Tacr3")

# now get all the probes and fold change data
# For the sake of space, can only show genes in the figure
# with fold change of greater than 2
fc.cutoff = 1
probes = unique(unlist(MsAnn$Symbol$Affy[genes]))
lfc = cyberT.results$logFC[probes,]
probes.2 = probes[ apply(abs(lfc)>log2(fc.cutoff),1,any) ]
lfc = cyberT.results$logFC[probes.2,]

# now collapse the expression data
lfc.collapsed = collapseExprMatrix( expMat=lfc,
                                    mapping=MsAnn$Symbol$Affy[genes],
                                    func=function(x) mean(abs(x)) )
lfc.symbols = lfc.collapsed[,c(1,3,2,4,6,5,7,8)]

# order the genes first by group. Then within each group, order them by 
# expression level at 16 hours
genes = rownames(lfc.symbols)
genegroups = numeric(length(genes))
for(i in 1:length(genelist)) { 
  genegroups = genegroups + i*(genes %in% genelist[[i]])
}
downgroups = 1-(!(genegroups %in% 9:10))*2
genes = genes[ order(genegroups,
                     lfc.symbols[,"A16-Sham16"]*downgroups ) ]


# Another type of display
fff =sign(lfc.symbols) * 2^abs(lfc.symbols)
aa = melt(fff,c("genes","treatment"))
library(stringr)
aa$toxin = str_match(aa$treatment, "[AB]*")
aa$time = as.numeric(str_match(aa$treatment,"[0-9]+"))
ggplot(data=aa, aes(x=time,y=value, group=toxin, color=toxin) ) + 
  scale_x_discrete(limits=colnames(fff),breaks=NULL) +
  geom_line() + facet_wrap( ~ genes ) + geom_point()

### Now the expression values
load("./data/RData/esets.RData")
eset = esets$all
ee = exprs(eset[probes,])
ee = collapseExprMatrix( expMat=ee,
                                   mapping=MsAnn$Symbol$Affy[genes],
                                   func=function(x) mean(abs(x)) )
pp = pData(eset)
mm = melt(ee, id.vars=pp$Toxin)
mm$Toxin = pp$Toxin[ match( mm$X2, rownames(pp) ) ]
mm$Time = pp$Time[ match( mm$X2, rownames(pp) ) ]
mm2 = melt(cast(mm,Toxin+Time~X1, mean))

ggplot(mm, aes(x=Time,y=value,group=interaction(Toxin),color=Toxin)) + 
         bgeom_jitter(position=position_jitter(width = .5)) +  
         geom_line(data=mm2,aes(x=Time,y=value,group=Toxin)) + 
         facet_wrap( ~X1, nrow=2 )
         


















