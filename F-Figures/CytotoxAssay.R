## prepare the data
library(signal)
library(reshape)
library(ggplot2)
library(plyr)

################################################################################
#                     From tables to objects                                   #
################################################################################

dat = read.csv("./data/CytotoxData.txt",sep="\t",skip=1)
conc = c(1000,100,1000,10,
  100,10,100,1,
  0,1,0,0.1,
  10,0.1,10,0.01,
  1,0.01,1,0.001,
  0.1,0.001,0.1,0.0001,
  0.01, 0.0001, 0.01, 0.00001,
  0.001, 0.00001, 0.001, 0 )
tox = c("A","B","A","B",
  "A","B","A","B",
  "C","B","C","B",
  "A","B","A","B",
  "A","B","A","B",
  "A","B","A","B",
  "A","B","A","B",
  "A","B","A","C" )
rows = rep( LETTERS[1:8],times=rep(4,8) )
cols = rep( 3:6, 8 )

dat$Time.Interval = NULL
colnames(dat) = c("Time",paste(tox,conc,rows,cols,sep="_") )

# the sweep numbers right before toxin was added for each column
firstsweep = 405
sweeps = c(46,75,53,83)
colTimes = dat$Time[firstsweep+sweeps]

# now melt the data to work with it easier
raw = melt(dat,id.vars="Time")
aa = strsplit( as.character(raw$variable), "_" )
bb = laply(aa,identity)
colnames(bb) = c("Toxin","Conc","Row","Col")
raw2 = data.frame( Time=raw$Time, bb, well=paste(bb[,"Row"],bb[,"Col"],sep=""),
                   value=raw$value )

# Correct the different times the toxin was added
map = setNames(colTimes, 3:6)
raw2$Time = raw2$Time - map[raw2$Col]
raw3 = raw2

# Interpolate all the wells to each other
reftime = raw3$Time[ raw3$well == "A3" ]
wells = levels(raw3$well)
for( well in wells ) {
  newdat = raw3[ raw3$well == well, c("Time","value") ]
  interpdat = interp1(newdat[[1]],newdat[[2]],reftime,"spline")
  raw3[ raw3$well == well, c("Time","value") ] = cbind(reftime,interpdat)
}

# Now normalize all of the wells to the zero timepoint
raw4 = raw3
norms = rep(raw4$value[ raw4$Time == 0 ], time=as.numeric(table(raw4$well)))
raw4$value = raw4$value/norms

# Now cast all of the data
foosd = function(x) ifelse(length(x)>1,sd(x),x*0)
means = cast(raw4, Time ~ Toxin+Conc, c(mean,foosd) )

# calculate the means and standard deviations
tp = melt(means,id.vars="Time")
tp2 = subset( tp, result_variable=="mean" )
tp2$sd = subset( tp, result_variable=="foosd" )$value

# finally plot the data
ggplot( tp2, aes(x=Time,y=value,group=interaction(Toxin,Conc)) ) +
  geom_ribbon( aes(ymin=value-sd,ymax=value+sd), fill=gray(0.5,0.3) ) + 
  geom_line( aes(colour=Conc) ) +
  xlim(c(-3,80)) + ylim(c(-0.1,1.8)) +
  facet_wrap( ~ Toxin )














