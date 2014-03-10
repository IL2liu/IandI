library(ggplot2)
library(reshape)
library(gridExtra)
library(plyr)
load("./data/annotation.RData/MsAnn.RData")


################################################################################
#                            Dose response                                     #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ raw$ID<200, c("ID","Toxin","ToxinDose",
             "Inflammation","MucosaThickness","CongestionExudate",
             "ArchitectureErosions","SubmucosalEdema","TotalPathology")]
mraw = melt(raw,id.vars=c("ID","Toxin","ToxinDose"))

# Make plot
ggplot( data=mraw, aes(x=factor(ToxinDose),y=value) ) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point( position=position_jitter(width=0.25,height=0) ) +
  facet_wrap( ~ variable, scale="free" )

################################################################################
#                     Cxcl1/2 mRNA expression data                             #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ !is.na(raw$Cxcl2.qRTPCR) &
            raw$Antibody!="None", 
           c("ID","Toxin","Time","Antibody","Death",
             "Cxcl1.qRTPCR","Cxcl2.qRTPCR" ) ]
raw$Death = factor(as.numeric(raw$Death!=""))
mraw = melt(raw,id.vars=c("ID","Toxin","Time","Antibody","Death"))

# Make plot
ggplot(mraw, aes(x=interaction(Antibody,Toxin), y=value) ) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point( position=position_jitter(width=0.25,height=0),
              aes(color=Death,shape=Death)) +
  facet_grid( variable ~ Time, scales="free" ) + 
  scale_color_manual(values=c("blue","red"))

################################################################################
#                   MPO staining after neutralization                          #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ !is.na(raw$MPOPositivesPerField) &
             raw$Antibody!="None", 
           c("ID","Toxin","Time","Antibody","Death",
             "MPOPositivesPerField" ) ]
raw$Death = factor(as.numeric(raw$Death!=""))


# Make plot
ggplot(raw, aes(x=interaction(Antibody,Toxin), y=MPOPositivesPerField) ) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point( position=position_jitter(width=0.25,height=0),
              aes(color=Death,shape=Death)) +
  scale_color_manual(values=c("blue","red"))

################################################################################
#                     Subcategories of pathology (time course)                 #
################################################################################

# the plots
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ raw$Antibody == "None" &
             raw$ID>199 &
             !is.na(raw$TotalPathology) , c("ID","Time","Toxin","TotalPathology",
                                            "Inflammation","MucosaThickness",
                                            "CongestionExudate","ArchitectureErosions",
                                            "SubmucosalEdema")]
mraw = melt(raw,id.vars=c("ID","Time","Toxin"))

ggplot( data=mraw, aes(x=interaction(Toxin,Time),y=value) ) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point( position=position_jitter(width=0.25,height=0) ) +
  facet_wrap( ~ variable, scale="free" )

### all pairwise p-values
# split data according to sample groups
groups = cast(mraw, ID~. | Time + variable + Toxin )

# reformat to contain only data values
for( i in 1:length(groups) ) {
  groups[[i]] = lapply( groups[[i]], 
                        function(x) lapply(x,function(y) y[,2]) )
}

# function to calculate pairwise p-values
pairwise.p = function(x) {
  pvals = matrix(1,length(x),length(x))
  dimnames(pvals) = list(names(x),names(x))
  for( i in 1:length(x) ) {
    for( j in 1:length(x) ) {
      if( is.null(x[[i]]) || is.null(x[[j]]) ) {
        pvals[i,j] = NA
      } else {
        pvals[i,j] = wilcox.test(x[[i]],x[[j]])$p.value
      }
    }
  }
  pvals[is.nan(pvals)] = 1
  pvals[ !apply(pvals,2,function(x) all(is.na(x))),
         !apply(pvals,1,function(x) all(is.na(x))) ]
}

# pairwise p-value for each measurement at each time
pvals = groups
for( i in 1:length(groups) ) {
  pvals[[i]] = lapply( groups[[i]], pairwise.p )
}

# only look at the Sham comparisons
sham.pvals = list()
for( i in 1:length(pvals) ) {
  sham.pvals[[i]] = sapply( pvals[[i]], function(x) x[rownames(x)!="Sham","Sham"] )
}
sham.pvals
lapply(sham.pvals,function(x) round(x,3))


################################################################################
#                           Blood counts from 300 and 900                      #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
raw = raw[ !is.na(raw$NE) &
             raw$Antibody=="None" &
             raw$ID > 199, 
           c("ID","Time","Toxin","Death",
             "WBC","NE","LY","MO","EO","BA","NE.p","LY.p","MO.p","EO.p",
             "Hb","HCT","MCV","MCH","MCHC","RDW","PLT","MPV") ]
raw$Death = factor(as.numeric(raw$Death!=""))
raw$Toxin = factor(raw$Toxin)
mraw = melt(raw,id.vars=c("ID","Time","Toxin","Death"))

# take out dead mice
mraw = mraw[ mraw$Death==0, ]
mraw$Death = NULL

# Plot all data
p = ggplot(data=mraw,aes(x=Toxin,y=value,fill=Toxin)) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point(position=position_jitter(height=0,width=0.25), aes(colour=Toxin)) +
  facet_wrap( ~ variable + Time, scale="free" ) + 
  scale_color_manual(values=c("blue","red","purple"))
p

# plot a selection of the data
vars = c("WBC","LY","NE" )
mraw2 = mraw[ mraw$variable %in% vars , ]
mraw2$variable = factor( droplevels(mraw2$variable), levels=vars )
p = ggplot(data=mraw2,aes(x=Toxin,y=value,fill=Toxin)) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point(position=position_jitter(height=0,width=0.25), aes(colour=Toxin)) +
  facet_grid( variable ~ Time, scale="free" ) + 
  scale_color_manual(values=c("blue","red","purple")) +
  theme( legend.position="none", axis.title=theme_blank() )
p1 = p

vars = c("MO","EO","BA")
mraw3 = mraw[ mraw$variable %in% vars, ]
mraw3$variable = factor( droplevels(mraw3$variable), levels=vars )
p2 = p %+% mraw3

vars = c("PLT","Hb","HCT" )
mraw3 = mraw[ mraw$variable %in% vars, ]
mraw3$variable = factor( droplevels(mraw3$variable), levels=vars )
p3 = p %+% mraw3

vars = c("MCV","MCH","MCHC","MPV","RDW","RBC")
mraw3 = mraw[ mraw$variable %in% vars, ]
mraw3$variable = factor( droplevels(mraw3$variable), levels=vars )
p4 = p %+% mraw3


grid.arrange(p1,p2,p3,ncol=3)


# do statistical tests on each comparison
allps = list()
alldata = cast(mraw, ID ~ . | Time + variable + Toxin, mean  )
for( time in names(alldata) ) {
  timedata = alldata[[time]]
  
  vars = names(timedata)
  varps = matrix(NA,nrow=length(vars),ncol=3)
  dimnames(varps) = list(vars,c("A.Sham","B.Sham","A.B"))
  for(var in names(timedata)) {
    vardata = timedata[[var]]
    vdata = lapply(vardata, function(x) x[[2]])
    toxins = names(vdata)
    
    fun = wilcox.test #wilcox.test
    p1 = fun( vdata$A, vdata$Sham )$p.value
    p2 = fun( vdata$B, vdata$Sham )$p.value
    p3 = fun( vdata$A, vdata$B )$p.value
    varps[var,] = c(p1,p2,p3)
  }
  allps[[time]] = varps
}
allps

# make it easier to see which groups meet significance cutoffs
sps = allps
for( i in 1:3 ) {
  
  ps = allps[[i]]
  tx = matrix("",ncol=ncol(ps),nrow=nrow(ps))
  dimnames(tx) = dimnames(ps)
  tx[ ps<0.01 ] = "**"
  tx[ ps<0.05 & ps >=0.01 ] = "*"
  tx[ ps<0.1 & ps >=0.05 ] = "-"
  sps[[i]] = data.frame(tx)
}
sps

################################################################################
#                    Colon pathology                                           #
################################################################################

# Parse data
raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
vars = c("Inflammation.Colon","MucosaThickness.Colon","CongestionExudate.Colon",
         "ArchitectureErosions.Colon","SubmucosalEdema.Colon","TotalPathology.Colon")
raw = raw[ !is.na(raw$NE) &
             raw$Antibody=="None" &
             raw$ID > 899, 
           c("ID","Time","Toxin",vars)]

mraw = melt(raw,id.vars=c("ID","Time","Toxin"))
ggplot( data=mraw, aes(x=interaction(Toxin,Time),y=value) ) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point( position=position_jitter(width=0.25,height=0) ) +
  facet_wrap( ~ variable, scale="free" )



















