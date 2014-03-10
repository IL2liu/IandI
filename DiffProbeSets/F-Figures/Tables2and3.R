load("./DiffProbeSets/data/RData/cameraMod.RData")
load("./DiffProbeSets/data/RData/cactus.RData")
library(GO.db)
library(qvalue)

################################################################################
#                       Helper functions                                       #
################################################################################

# A versioin of the qvalue function that handles NA values by returning
# NA values
qvalue.na = function(pvals) {
  out = pvals*0
  qvals = qvalue(na.omit(pvals))$qvalue
  out[names(qvals)] = qvals
  return(out)
}

################################################################################
#                               Table 2                                        #
################################################################################

# only use some of the gene set groups
cac.result = cactus.result[c("mf","bp","cc")]

table2 = list()
for( comparison in c("A2-Sham2","B2-Sham2") ) {
# Extract p-values for comparisons at 2hrs (for multiple gene set groups)
  pvals = lapply(cac.result, function(x) x[,comparison])
  
  # get q-values from these p-values
  qvals = lapply(pvals, function(x) qvalue.na(x) )
  
  # put everything into a matrix that is easier to handle
  db.names = rep(names(qvals),times=sapply(qvals,length))
  set.names = as.character(unlist( lapply(qvals,names) ))
  all.qvals = as.numeric(unlist(qvals))
  all.pvals = as.numeric(unlist(pvals))
  scores = data.frame(set=set.names, pval=all.pvals, qval=all.qvals, db=db.names )
  
  # get which sets meet a qvalue<0.1 cut off and sort
  sigs = na.omit( scores[ scores$qval<0.05, ] )
  
  # now just some formatting
  sigs$set = apply(sigs, 1, function(x) ifelse( x["db"] %in% c("mf","bp","cc"),
                                                Term(x["set"]), x ) )
  sigs$logpval = -log10(sigs$pval)
  table2[[comparison]] = sigs[ order(sigs$db,-sigs$logpval), c("set","logpval","qval","db")]
}
table2

################################################################################
#                               Table 3                                        #
################################################################################

# The top gene sets for A16 versus Sham16
pvals = camera.result$mf[,"A16-Sham16"]
qvals = qvalue(pvals)$qvalues

set.qvals = sort( qvals[ qvals<0.3 ] )
set.names = Term(names(set.qvals))
set.pvals = pvals[names(set.qvals)]


tab = data.frame( name=substr(set.names,1,50), 
                  logpval=-log10(set.pvals), 
                  qval=set.qvals )
table3 = tab[ order(set.pvals), ]
table3




















