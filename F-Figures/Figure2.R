library(ggplot2)
library(reshape)


################################################################################
#                       Histopathology of all mice                             #
################################################################################

raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ raw$Antibody == "None" &
           raw$ID>199 &
           !is.na(raw$TotalPathology) , c("ID","Time","Toxin","TotalPathology")]

ggplot(data=raw,aes(x=interaction(Toxin,Time),y=TotalPathology,fill=Toxin)) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point(position=position_jitter(height=0,width=0.25))

# Statistical tests -- nonparametric
dat = cast(raw,ID~Toxin+Time,value="TotalPathology")
dat = lapply(dat,function(x) x[!is.na(x)])

# at 2 hours
wilcox.test(dat$A_2,dat$Sham_2) # 0.003318  *
wilcox.test(dat$AB_2,dat$Sham_2) # 0.02144  -
wilcox.test(dat$B_2,dat$Sham_2) # 0.2194    -
wilcox.test(dat$A_2,dat$B_2) # 0.01276   *
wilcox.test(dat$A_2,dat$AB_2) # 0.5169   -

# at 6 hours
wilcox.test(dat$A_6,dat$Sham_6) # 0.004922   *
wilcox.test(dat$AB_6,dat$Sham_6) # 0.896    -
wilcox.test(dat$B_6,dat$Sham_6) # 0.6847     -
wilcox.test(dat$A_6,dat$B_6) # 0.007796   *
wilcox.test(dat$A_6,dat$AB_6) # 0.05182   *

# at 16 hours
wilcox.test(dat$A_16,dat$Sham_16) # 0.000118  *
wilcox.test(dat$B_16,dat$Sham_16) # 0.03227  *
wilcox.test(dat$A_16,dat$B_16) # 0.001577   *


################################################################################
#                     MPO for experiments 200 and 300                          #
################################################################################

raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ raw$Antibody == "None" &
             raw$ID>199 &
             !is.na(raw$MPOPositivesPerField) , c("ID","Time","Toxin","MPOPositivesPerField")]

ggplot(data=raw,aes(x=interaction(Toxin,Time),y=MPOPositivesPerField,fill=Toxin)) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point(position=position_jitter(height=0,width=0.25))


# Statistical tests -- nonparametric
dat = cast(raw,ID~Toxin+Time,value="MPOPositivesPerField")
dat = lapply(dat,function(x) x[!is.na(x)])

# at 2 hours
wilcox.test(dat$A_2,dat$Sham_2) # 0.7309
wilcox.test(dat$AB_2,dat$Sham_2) # 0.1858
wilcox.test(dat$B_2,dat$Sham_2) # 1
wilcox.test(dat$A_2,dat$B_2) # 0.164
wilcox.test(dat$A_2,dat$AB_2) # 0.2

# at 6 hours
wilcox.test(dat$A_6,dat$Sham_6) # 0.1255
t.test(dat$A_6,dat$Sham_6) # 0.05548
wilcox.test(dat$AB_6,dat$Sham_6) # 1
wilcox.test(dat$B_6,dat$Sham_6) # 0.7922
wilcox.test(dat$A_6,dat$B_6) # 0.1797
wilcox.test(dat$A_6,dat$AB_6) # 0.381

# at 16 hours
wilcox.test(dat$A_16,dat$Sham_16) # 0.002953
wilcox.test(dat$B_16,dat$Sham_16) # 0.02
wilcox.test(dat$A_16,dat$B16) # 0.007813

################################################################################
#                           Venn diagrams                                      #
################################################################################

load("./data/RData/cyberT.results.RData")
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

load("./data/RData/cyberT.results.RData")
library(gplots)

logFC = cyberT.results$logFC
colnames(logFC) = c("A2","B2","AB2","A6","B6","AB6","A16","B16")
logFC = logFC[,c("B2","A2","AB2","B6","A6","AB6","B16","A16")]
cors = cor(logFC)

heatmap.2(cors,col=bluered,symbreak=TRUE,Rowv=NULL,
          Colv=NULL,dendrogram="none",trace="none",
          cellnote=round(cors,2),breaks=1000,
          density.info="none",asp=1)










