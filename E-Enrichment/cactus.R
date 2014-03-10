source("./C-Annotations/CollapseIDs.R")
source("./D-DiffExpression/cybertWrapper.R")
source("./E-Enrichment/indexGeneSets.R")
require(limma)

# Null hypothesis. Each gene's test statistic has a standard normal. Hence,
# the expected sum is 0, and the expected variance is 
# the sum of all the elements in the covariance matrix (i.e. n + n(n-1)p where
# p is the average correlation and n is the gene set size)
cactus = function( eset, genesets, comparisons, mapping=NULL ) {

  # get sample groups
  data = exprs(eset)
  pdat = pData(eset)
  groups = interaction( pdat$Toxin, pdat$Time, sep="" )
  
  # prepare matrix for the output
  out = matrix(NA,length(genesets),length(comparisons))
  dimnames(out) = list( names(genesets), 
                          lapply(comparisons, function(x) paste(x,collapse="-")) )
  
  for( i in seq_along(comparisons) ) {
    
    # separate sample groups
    comparison = comparisons[[i]]
    groupE.label = comparison[1]
    groupC.label = comparison[2]
    data.E = as.matrix( data[ , groups == groupE.label ] )
    data.C = as.matrix( data[ , groups == groupC.label ] )
    data.combined = cbind(data.E,data.C)
    
    # collapse the data matrix if a mapping is provided
    if( !is.null(mapping) ) {
      
      result = collapseExprMatrix( data.combined, mapping, return.probes=TRUE)
      new.eset = new("ExpressionSet",exprs=result$expMat)
      
      # and correct the gene sets after collapsing the matrix
      indexed.genesets = indexGeneSets( genesets, result$collapseMap, 10, 700 )
    } else {
      
      new.eset = new("ExpressionSet",exprs=data.combined)
      names(genes) = genes = rownames(exprs(new.eset))
      indexed.genesets = indexGeneSets( genesets, genes, 10, 700 )
    }
    pData(new.eset) = pData(eset)[ match( colnames(exprs(new.eset)), 
                                          rownames(pData(eset)) ), ]
    
    # get the t-statistics
    conf=10
    cybert = cybertWrapper(new.eset,list(comparison),conf)
    tstats = cybert$tstats[,1]
    deg.free = cybert$degfree[1,1] # bayes/posterior degrees of freedom
    
    #### This starts the actual code for the enrichment ####
    
    # get z scores in the same quantile as the t distribution
    zstats = zscoreT(tstats,deg.free)
    
    # get the average intergene correlation for each gene set
    aveCor = function(x) {
      coeff = cor(t(exprs(new.eset)[x,]))
      aver.cor = mean( coeff[ upper.tri(coeff) ] )
      return(aver.cor)
    }
    geneset.corr = sapply( indexed.genesets, aveCor )
    
    # the variance for each gene set adjusted for the correlation
    m = sapply(indexed.genesets,length)
    geneset.vars = (1/m) * ( 1 + (m-1)*geneset.corr )
    
    # the mean z statistic for each gene set
    geneset.means = sapply(indexed.genesets, function(x) mean(zstats[x]))
        
    pvals = geneset.means*0+1
    for( j in 1:length(indexed.genesets) ) {
      mu = geneset.means[j]
      s2 = geneset.vars[j]
      low = pnorm(mu, mean=0, sd=sqrt(s2))
      high = pnorm(-mu, mean=0, sd=sqrt(s2))
      two.sided = min(low,high)*2
      pvals[j] = two.sided
    }
    
    out[names(pvals),i] = pvals
  }
  
  return( out[apply(out,1,function(x) !all(is.na(x))),] )

}