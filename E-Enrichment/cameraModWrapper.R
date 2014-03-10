source("./C-Annotations/CollapseIDs.R")
source("./D-DiffExpression/cybertWrapper.R")
source("./E-Enrichment/indexGeneSets.R")
source("./E-Enrichment/cameraMod.R")

cameraModWrapper = function( eset, genesets, comparisons, mapping=NULL  ) {

  # prepare data matrix for output
  pvals = matrix(NA,length(genesets),length(comparisons))
  dimnames(pvals) = list( names(genesets), 
                          lapply(comparisons, function(x) paste(x,collapse="-")) )
  
  # get sample groups
  data = exprs(eset)
  pdat = pData(eset)
  groups = interaction( pdat$Toxin, pdat$Time, sep="" )
  
  # do an enrichment for each comparison
  for( i in seq_along(comparisons) ) {
    
    comparison = comparisons[[i]]
    
    # set up the data matrix
    groupE.label = comparison[1]
    groupC.label = comparison[2]
    data.E = as.matrix( data[ , groups == groupE.label ] )
    data.C = as.matrix( data[ , groups == groupC.label ] )
    data.combined = cbind(data.E,data.C)
    
    # set up the design matrix
    column1 = c(rep(1,ncol(data.E)),rep(0,ncol(data.C)))
    design = cbind( Intercept=rep(1,length(column1)), Control=1-column1 )
    colnames(design) = c(groupE.label, groupC.label)
    
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
    
    # run camera
    result = cameraMod(indices=indexed.genesets, exprs(new.eset), design, 
                       contrast=ncol(design),tstats,d0=2*conf)
    pvals[rownames(result),i] = result[,"TwoSided"]
  }
  return( pvals[apply(pvals,1,function(x) !all(is.na(x))),] )
}









