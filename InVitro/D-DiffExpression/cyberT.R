# MUST RUN IBMT.R first
source("./InVitro/D-DiffExpression/cyberTFunction.R")

# the expression set data
load("./InVitro/data/RData/bppsl.eset.RData")
eset = bppsl.eset

# cheat and get the matrices already made for IBMT
load("./InVitro/data/RData/ibmtResults.RData")
dd = ibmt.results$design
cm = ibmt.results$contrasts


################################################################################
#                  cyberT on batch corrected data set                          #
################################################################################

cybert.hs.results = cybert.wrapper( eset, dd, cm )
save(cybert.hs.results, file="./InVitro/data/RData/cybert.hs.results.RData")
toptable.ct( cybert.hs.results, 2, 20 )


################################################################################
#             cyberT on batch corrected, collapsed data set                    #
################################################################################

# steal collapsed data set from IBMT code
load("./InVitro/data/RData/collapsed.ibmtResults.RData")
eset = collapsed.ibmt.results$collapsed.eset
cybert.hs.collapsed.results = cybert.wrapper( eset, dd, cm )
save(cybert.hs.collapsed.results, file="./InVitro/data/RData/cybert.hs.collapsed.results.RData")
toptable.ct( cybert.hs.collapsed.results, 1, 30 )

################################################################################
#          cyberT on batch corrected, collapsed data sets                      #
################################################################################

# steal the collapsed data sets from the IBMT code
load("./InVitro/data/RData/mc.ibmt.RData")
mc.cybert = vector("list",ncol(cm))

# run IBMT for each collapsed contrast
for( i in 1:ncol(cm) ) {
  subeset = mc.ibmt[[i]]$collapsed.eset
  
  arrays = apply( dd[ , cm[,i]!=0 ]!=0, 1, any )
  subdd = dd[arrays,cm[,i]!=0]
  subcm = t(t(cm[cm[,i]!=0,i]))
  
  mc.cybert[[i]] = cybert.wrapper( subeset, subdd, subcm )
}

# now combine the multiple lists
combined = list()
for( i in names(mc.cybert[[1]]) ) {
  combined[[i]] = sapply( mc.cybert, function(x) x[[i]] )
  rownames(combined[[i]]) = rownames(mc.cybert[[1]][[1]])
  colnames(combined[[i]]) = colnames(cm)
}
mc.cybert = combined

save(mc.cybert,file="./InVitro/data/RData/mc.cybert.RData")











