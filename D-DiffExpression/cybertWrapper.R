# A very specific function. Runs cyber T for the comparisons given
# The labels for the comparison must be in pData(eset)
source("./D-DiffExpression/bayesreg.R")
library(affy)

cybertWrapper = function( eset, comparisons, conf=10 ) {
  
  # get sample groups
  pdat = pData(eset)
  groups = interaction( pdat$Toxin, pdat$Time, sep="" )
  
  # set up the matrices for the data and the results
  data = exprs(eset)
  pval = matrix(NA,nrow(data),length(comparisons))
  colnames(pval) = sapply(comparisons,function(x) paste(x,collapse="-"))
  rownames(pval) = rownames(data)
  sigC = sigE = sdC = sdE = bayesSDC = bayesSDE = rasdC = rasdE = degfree = tstats = logFC = pval
  
  # running cyberT for each comparison. Think of "E" as referring to
  # experimental group and "C" referring to control group
  for( i in seq_along(comparisons)) {
    comparison = comparisons[[i]]
    groupE.label = comparison[1]
    groupC.label = comparison[2]
    
    data.E = as.matrix( data[ , groups == groupE.label ] )
    data.C = as.matrix( data[ , groups == groupC.label ] )
    data.combined = cbind(data.C,data.E)
    
    result = bayesT(data.combined,numC=ncol(data.C),numE=ncol(data.E),ppde=FALSE,
                    betaFit=2, bayes=TRUE, winSize=101, conf=conf)
    
    pval[,i] = result[,"pVal"]
    logFC[,i] = result[,"fold"]
    tstats[,i] = result[,"bayesT"]
    degfree[,i] = result[,"bayesDF"]
    rasdC[,i] = result[,"rasdC"]
    rasdE[,i] = result[,"rasdE"]
    bayesSDC[,i] = result[,"bayesSDC"]
    bayesSDE[,i] = result[,"bayesSDE"]
    sigC[,i] = result[,"meanC"]
    sigE[,i] = result[,"meanE"]
    sdC[,i] = result[,"stdC"]
    sdE[,i] = result[,"stdE"]
  }
  
  return( list(pval=pval,logFC=logFC,
               tstats=tstats,degfree=degfree,
               rasdC=rasdC, rasdE=rasdE,
               bayesSDC=bayesSDC, bayesSDE=bayesSDE,
               sigC=sigC, sigE=sigE, sdC=sdC, sdE=sdE ) )
}





