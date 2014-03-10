source("./InVitro/D-DiffExpression/cyberTFunction.R")
library(affy)

# the expression set data
load("./InVitro/DiffProbeSets/data/RData/bppsl.eset.RData")
eset = bppsl.eset

# cheat and get the matrices already made for IBMT
load("./InVitro/DiffProbeSets/data/RData/ibmtResults.RData")
dd = ibmt.results$design
cm = ibmt.results$contrasts


################################################################################
#                  cyberT on batch corrected data set                          #
################################################################################

cybert.hs.results = cybert.wrapper( eset, dd, cm )
save(cybert.hs.results, file="./InVitro/DiffProbeSets/data/RData/cybert.hs.results.RData")
toptable.ct( cybert.hs.results, 1, 20, probeann="ENSG" )

################################################################################
#     cyberT on batch corrected, collapsed to Gene Symbol data set             #
################################################################################

# steal collapsed data set from IBMT code
load("./InVitro/DiffProbeSets/data/RData/collapsed.ibmtResults.RData") 
eset = collapsed.ibmt.results$collapsed.eset
cybert.hs.collapsed.results = cybert.wrapper( eset, dd, cm )
save(cybert.hs.collapsed.results, file="./InVitro/DiffProbeSets/data/RData/cybert.hs.collapsed.results.RData")
toptable.ct( cybert.hs.collapsed.results, 4, 30 )


