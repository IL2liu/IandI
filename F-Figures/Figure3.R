source("./F-Figures/genesWithin.R") # gtpasegenes, metgenes
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


################################################################################
#                            Figures 3A and 3B                                 #
################################################################################

## The list of the genes to be shown in the figure
genelist = list(ifn=gtpasegenes$ifn,
                small=gtpasegenes$small,
                gap=gtpasegenes$gap,
                gef=gtpasegenes$gef,
                G=gtpasegenes$G,
                tub=gtpasegenes$tub,
                rho.bind=gtpasegenes$rho.binding,
                uptake=metgenes$uptake,
                chol=metgenes$chol,
                FAmet=metgenes$FAmet , # down
                enzym=metgenes$enzymes ) # down
genes = unlist(genelist)

# now get all the probes and fold change data
# For the sake of space, can only show genes in the figure
# with fold change of greater than 2
fc.cutoff = 2
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

# prepare the data for ggplot 2
# in the svg conversion, the first and last rows' bubbles end up as polygons?
# so... making a dummy first and last rows
lfc.symbols = rbind(blank=rep(log2(61),8),
                    blank2=rep(log2(61),8),
                    lfc.symbols,
                    blank3=rep(log2(61),8),
                    blank4=rep(log2(61),8)) 
bb = melt(lfc.symbols,c("genes","treatment"))

## Now do the bubble charts
p = ggplot(data=bb, aes(x=treatment,y=genes) ) +
  scale_y_discrete(limits=c("blank3","blank4",rev(genes),"blank","blank2")) +
  scale_x_discrete(limits=colnames(lfc.symbols),breaks=NULL) +
  theme_bw() +
  theme(panel.grid.major=element_line(color="gray"),
        axis.ticks=element_blank(),
        panel.border=element_blank(),
        legend.position="right",
        axis.title=element_blank(),
        axis.text.y=element_text(size=6.5) )

#devSVG(file="./gtpases.svg",height=23,width=6)
p + geom_point(aes(size=abs(value),color=(value>0)*2-1)) +
    scale_size_area(breaks=log2(c(1,1.2,2,5,15,50)),
                    limits=c(0,max(abs(lfc.symbols))*1.05)) +
    scale_color_gradient(low="blue",high="red")
#dev.off()
# format the svg with Illustrator

################################################################################
#                              Figure 3C                                       #
################################################################################

# Gene ontology inflammation associate
gos = c("GO:0006955","GO:0045087","GO:0045088","GO:0045089","GO:0006954",
        "GO:0050729","GO:0006935","GO:0030593","GO:0008009","GO:0005125",
        "GO:0019221","GO:0004896","GO:0019955","GO:0050710","GO:0035458",
        "GO:0034341","GO:0006959","GO:0070498","GO:0005149")
db = c(MsAnnGO$BP$MGI,MsAnnGO$MF$MGI,MsAnnGO$CC$MGI)

# Get the genes in these go categroies
mgis = unique(as.vector(unlist(db[gos])))
m2p = MsAnn$MGI$Affy[mgis]
m2p = m2p[sapply(m2p,length)>0]
goprobes = unique(as.vector(unlist(m2p)))
gomgis = mgis

# REACTOME inflammation associate
reactsets = c( "REACTOME_INNATE_IMMUNITY_SIGNALING",
               "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
               "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
               "REACTOME_SIGNALING_IN_IMMUNE_SYSTEM")
entrezs = unique(unlist(GSAnn$REACTOME$HumanENTREZ[reactsets]))
rprobes = unique(unlist(MsAnn$HumanENTREZ$Affy[entrezs]))
rmgis = unique(unlist(MsAnn$Affy$MGI[rprobes]))

probes = union(goprobes,rprobes)
mgis = union(gomgis,rmgis)

# Use MGI symbols to collapse the probe sets to genes
summary( probes %in% rownames(cyberT.results$logFC) )
probes2 = intersect( probes, rownames(cyberT.results$logFC) )
hm = cyberT.results$logFC[ probes2, ]
hm.mgi = collapseExprMatrix(hm,MsAnn$MGI$Affy)
rownames(hm.mgi) = sapply( MsAnn$MGI$Symbol[rownames(hm.mgi)], function(x) x[1] )
hm2 = hm.mgi

#################### Now choose the clusters of genes

# the clusters of genes
g1 = c("Cxcl1","Cxcl2","Cxcl3","Cxcl10","Ccl3","Ptgs2") # red
g2 = c("S100a9","Reg3g","S100a8","Areg","Chi3l3","Reg3b","Il1b","Tmem173","Tac1",
       "Fpr2","Il1rn","Adam8","Tnfaip3","Srgn","Tnfsf9","Cd55","Clec4d",
       "Chst4","Dusp6","Icam1","Tnf") # blue
g3 = c("Ido1","Il33","Iigp1","Gbp2","Lcn2","Ereg","Igtp","Ecscr","Cd177","Ecm1",
       "Tgtp1","Cxcl9","Gbp6") # green
g4 = c("Nlrc5","Slc7a9","Cd274","Tgm2","Ifitm1","C3","Irgm2","Irgm1","Reg3a",
       "Zbp1","Gbp3","Gata3","Serping1","Gm4951","Gbp10","Osmr") # gold
g5 = c("Edn2","Bmp2","Mertk","Grb7","Il18","Gas6","Ephx2","Il1r1","Bmp8a",
       "Edn1","Il15","Cmtm8","Prlr","Agxt2","Afap1l2","Map2k6","Pycard") # tarheel blue
duplicate.genes = c("Gm6522","Tgtp2","Iigp1b","Gbp7")
gother = setdiff( rownames(hm2), c(g1,g2,g3,g4,g5,duplicate.genes) )
genelist = list(g1,g2,g3,g4,g5,gother)
colors = rep(1:length(genelist),times=sapply(genelist,length))
hm3 = data.frame(hm2[ unlist(genelist), ])
hm3$groups = factor(colors)
hm3 = hm3[rev(rownames(hm3)),]

##################### Now color the groups of genes and plot

# Toxin A 6h versus 16h
p = ggplot(hm3,aes(x=A6.Sham6,y=A16.Sham16)) + 
  geom_abline(intercept=0,slope=1,color="gray") +
  geom_point(alpha=0.8,aes(color=groups)) +
  geom_text(label=rownames(hm3),size=1) +
  scale_y_continuous(breaks=2*(-1:3),limits=c(-3.5,6.8)) + 
  scale_x_continuous(breaks=c(-2,0,2,4),limits=c(-2.3,4.9)) +
  coord_fixed() +
  scale_color_manual(values=c("red","blue","green","orange","lightblue","black"))

# Toxin B 6h versus 16h
q = p
q$mapping = aes(x=B6.Sham6,y=B16.Sham16)
grid.arrange(p,q,nrow=1)


################ Now do line plots for each cluster

# Toxin A plots of each group
hm4 = hm3[ hm3$groups!=6, ]
hm4$groups = factor(hm4$groups)
lfc = hm4[,1:8 ]
fc = sign(lfc)*(2^abs(lfc))
fc.zeroed = fc - sign(fc)
hm4[,1:8] = fc.zeroed
hm4$AB0 = hm4$B0 = hm4$A0 = 0
hm4$gene = rownames(hm4)
mraw = melt(hm4,c("groups","gene"))

timing = toxin = mraw$variable
levels(timing) = c(2,2,2,6,6,6,16,16,0,0,0)
levels(toxin) = c("A","B","AB","A","B","AB","A","B","A","B","AB")
mraw$time = as.numeric(as.character(timing))
mraw$toxin = toxin

araw = mraw[ mraw$toxin %in% c("A","B"),  ]
r = ggplot(araw,aes(x=time,y=value,group=gene)) +
    geom_line(size=0.25*.pt) + 
    geom_point(size=0.8) +
    geom_text(data=araw[araw$time==16,],aes(x=16.3,y=value,label=gene),
              size=5*.ptText,hjust=0) +
    facet_wrap( groups ~ toxin , scales="free_y",ncol=2)
r

# Combined chemokine plot with toxin a and toxiin b
r2 = ggplot(araw,aes(x=time,y=value,group=interaction(gene,toxin),color=toxin)) +
  geom_line(size=0.25*.pt) + 
  geom_point(size=0.8) +
  geom_text(data=araw[araw$time==16,],aes(x=16.3,y=value,label=gene),
            size=5*.ptText,hjust=0) +
  facet_wrap( ~ groups , scales="free_y",ncol=2)
r2

################################################################################
#                 For statements in results and discussion                     #
################################################################################

# Statements about correlation and magnitude of expression change
hm4 = hm3[ hm3$groups %in% 2:5, ]
cors = cor(hm4[,-ncol(hm4)])
cors["A2.Sham2","B2.Sham2"]
cors["A6.Sham6","B6.Sham6"]
cors["A16.Sham16","B16.Sham16"]

hm5 = hm3[ hm3$groups %in% 1:5, ]
f2 = as.matrix(hm5[,"A2.Sham2"] - hm5["B2.Sham2"])
f6 = as.matrix(hm5[,"A6.Sham6"] - hm5["B6.Sham6"])
f16 = as.matrix(hm5[,"A16.Sham16"] - hm5["B16.Sham16"])

uf2 = sign(f2)*2^abs(f2) - sign(f2) # unlog the fold change
mean(uf2) # 0.22
uf6 = sign(f6)*2^abs(f6) - sign(f6) # unlog the fold change
mean(uf6) # 2.32
uf16 = sign(f16)*2^abs(f16) - sign(f16) # unlog the fold change
mean(uf16) # 3.05







