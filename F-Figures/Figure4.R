library(ggplot2)
library(reshape)

################################################################################
#                              Figure 4A - Protein                             #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ !is.na(raw$Cxcl1.Serum) &
           !is.na(raw$Cxcl2.Serum) , 
           c("ID","Time","Toxin","Antibody","Death","Cxcl1.Serum","Cxcl2.Serum")]
raw$Death = factor(as.numeric(raw$Death!=""))
mraw = melt(raw,measure.vars=c("Cxcl1.Serum","Cxcl2.Serum"))

# Plot
p = ggplot(data=mraw,aes(x=interaction(Antibody,Toxin),y=value,fill=Toxin)) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point(position=position_jitter(height=0,width=0.25),
             aes(color=Death,shape=Death)) +
  facet_wrap( ~ variable ) + 
  scale_color_manual(values=c("blue","red"))
p

# Statistical tests -- nonparametric
dat = cast(mraw,ID~Toxin+Antibody+variable,value="TotalPathology")
dat = lapply(dat,function(x) x[!is.na(x)])
wilcox.test( dat[["A_Anti-Cxcl_Cxcl1.Serum"]], dat$A_Isotype_Cxcl1.Serum ) # p=0.00236
wilcox.test( dat[["A_Anti-Cxcl_Cxcl2.Serum"]], dat$A_Isotype_Cxcl2.Serum ) # p=0.01998


################################################################################
#                     Figure 4B - Neutrophil count                             #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ !is.na(raw$NE) &
            raw$Antibody!="None", 
           c("ID","Time","Toxin","Antibody","Death","NE")]
raw$Death = factor(as.numeric(raw$Death!=""))
mraw = melt(raw,measure.vars=c("NE"))

# Plot
p %+% mraw + facet_wrap( ~ Time )

# Statistical tests -- nonparametric
dat = cast(mraw,ID~Toxin+Antibody+Time,value="NE")
dat = lapply(dat,function(x) x[!is.na(x)])
wilcox.test( dat[["Sham_Isotype_6"]], dat$A_Isotype_6 )
wilcox.test( dat[["A_Anti-Cxcl_6"]], dat$A_Isotype_6 ) # p=0.3823
wilcox.test( dat[["A_Anti-Cxcl_16"]], dat$A_Isotype_16 ) # p=0.05714
2/choose(7,3) # minimum p-value=0.05714
t.test( dat[["A_Anti-Cxcl_6"]], dat$A_Isotype_6 ) # p=0.603
t.test( dat[["A_Anti-Cxcl_16"]], dat$A_Isotype_16 ) # p=0.01308


################################################################################
#                        Figure 4C - Pathology                                 #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ !is.na(raw$TotalPathology) &
             raw$Antibody!="None", 
           c("ID","Time","Toxin","Antibody","Death","TotalPathology")]
raw$Death = factor(as.numeric(raw$Death!=""))
mraw = melt(raw,measure.vars=c("TotalPathology"))

# Plot
p %+% mraw + facet_wrap(~Time)

# Statistical tests -- nonparametric
dat = cast(mraw,ID~Toxin+Antibody+Time,value="NE")
dat = lapply(dat,function(x) x[!is.na(x)])
wilcox.test( dat[["A_Anti-Cxcl_6"]], dat$A_Isotype_6 ) # p=0.5624
wilcox.test( dat[["A_Anti-Cxcl_16"]], dat$A_Isotype_16 ) # p=1


################################################################################
#                          Figure 4D - Survival                                #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ raw$Antibody!="None", 
           c("ID","Time","Toxin","Antibody","Death")]
raw$Death = factor(as.numeric(raw$Death!=""))

# Survival
s = cast( raw, Death~Time+Toxin+Antibody, value="ID", fun=length )
apply(s,2,function(x) paste(x[1],"/",sum(x)))

















