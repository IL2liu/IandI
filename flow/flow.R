library(reshape)
library(ggplot2)
library(gridExtra)

aa = read.csv("./data/flowData/flowData.csv",sep="\t")
bb = aa

# Panel A
# Percentage of epithelial and non-epithelial
cc = bb[ bb$Panel=="A", ]
cc$Epithelial = 100-cc$Nonepithelial/cc$Singlets*100
cc$CD3 = cc$CD3 / cc$CD45 * 100
cc$CD45 = cc$CD45 / cc$Nonepithelial * 100
cc$Nonepithelial = cc$Nonepithelial / cc$Singlets * 100

dd = cc[,c(1,5,10,6,7,7)]
dd$Toxin = c("Sham","Sham","Sham","B","B","B","A","A")
araw = melt(dd,id.vars=c("ID","Toxin"))
araw$ID = factor(araw$ID)

mytheme = theme( panel.grid.minor = element_blank(),
                 panel.grid.major.x = element_blank(),
                 legend.position="none",
                 strip.background=element_blank(),
                 axis.ticks=element_blank(),
                 axis.title=element_blank(),
                 axis.text.x=element_blank())
a = ggplot( data=araw, aes(x=Toxin,y=value) ) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point( aes(shape=ID,color=Toxin), size=3.2 ) +
  facet_wrap( ~ variable, scale="free" ) +
  scale_shape_manual(values=c(16:18,16:18,16:17)) +
  scale_x_discrete(limits=c("Sham","A","B")) + 
  mytheme

# Do t-tests to see if there are differences
pp = list()
for( type in c("Epithelial","Nonepithelial","CD45","CD3") ) {
  bc = t.test(dd[1:3,type],dd[4:6,type])
  ac = t.test(dd[1:3,type],dd[7:8,type])
  ab = t.test(dd[4:6,type],dd[7:8,type])
  pp[[type]] = c( ac=ac$p.value, bc=bc$p.value, ab=ab$p.value )
}
pp1=pp

# Panel B
cc = bb[ bb$Panel=="B", ]
cc$CD11b = cc$CD11b / cc$CD45 * 100
cc$B220 = cc$B220 / cc$CD45 * 100
cc$CD45 = cc$CD45 / cc$Nonepithelial * 100
cc$Epithelial = 100-cc$Nonepithelial/cc$Singlets*100
cc$Nonepithelial = cc$Nonepithelial / cc$Singlets

dd = cc[, c(1,5,10,6,8,9)]
dd$Toxin = c("Sham","Sham","Sham","B","B","B","A","A")
braw = melt(dd,id.vars=c("ID","Toxin"))
braw$ID = factor(braw$ID)

b = ggplot( data=braw, aes(x=Toxin,y=value) ) +
  stat_summary(fun.y="mean",geom="bar") +
  geom_point( aes(shape=ID,color=Toxin), size=3.2 ) +
  facet_wrap( ~ variable, scale="free" ) +
  scale_shape_manual(values=c(16:18,16:18,16:17)) +
  scale_x_discrete(limits=c("Sham","A","B")) +
  mytheme

pp = list()
for( type in c("Epithelial","Nonepithelial","CD45","CD11b","B220") ) {
  bc = t.test(dd[1:3,type],dd[4:6,type])
  ac = t.test(dd[1:3,type],dd[7:8,type])
  ab = t.test(dd[4:6,type],dd[7:8,type])
  pp[[type]] = c( ac=ac$p.value, bc=bc$p.value, ab=ab$p.value )
}
pp2 = pp

grid.arrange(a,b)


