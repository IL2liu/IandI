source("./D-DiffExpression/bayesreg.R")

cybert.wrapper = function( eset, dd, cm ) {
  
  # run the test on each contrast
  cts = vector("list",ncol(cm))
  for( i in 1:ncol(cm) ) {
    exp.arrays = dd[ , cm[,i]!=0 ][,1] != 0
    con.arrays = dd[ , cm[,i]!=0 ][,2] != 0
    data.com = cbind(exprs(eset[,con.arrays]),exprs(eset[,exp.arrays]))
    cts[[i]] = bayesT( data.com, sum(con.arrays), sum(exp.arrays), ppde=FALSE, betaFit=1 )
  }
  
  # format the output
  out=list()
  for( i in c("fold","bayesT","bayesDF","pVal") ) {
    out[[i]] = sapply(cts,function(x) x[,i])
    colnames(out[[i]]) = colnames(cm)
  }
  
  return(out)
}

toptable.ct = function( ctres, i=1, n=20, probeann="Affy" ) {
  if(!exists("HsAnn")) {load("./InVitro/data/annotation.RData/HumanAnnotations.RData")}
  pv = ctres$pVal[,i]
  probesets = names(sort(pv)[1:n])
  lfc = ctres$fold[probesets,i]
  if(probeann=="Affy") {
    syms=sapply( HsAnn$Affy$Symbol[probesets], function(x) paste(x,collapse=":::"))
  } else if(probeann=="ENSG") {
    syms = sapply( HsAnn$Ensembl$Symbol[probesets], function(x) paste(x,collapse=":::") )
  }  
  pvs = pv[probesets]
  tab = data.frame(syms,lfc,pvs)
  return(tab)
}

