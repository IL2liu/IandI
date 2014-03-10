# Very often, Affymetrix probe sets map to multiple gene-level identifiers (e.g. 
# ENTREZ, Symbol, Ensembl). Their are multiple ways to collapse the
# multi-mapping probe sets into one gene. With a single vector of values (e.g. 
# p-values or a test statistic), estimators such as the max, mean, median may 
# be used. With a matrix of values, combining the probes becomes more complicated.
# For instance, if one takes the probe set with the max intensity in sample 1 (column 1)
# of the matrix, then what do they do if another probe set has the maximum
# intensity for sample 2? Hence, one may choose to use only one of the probe sets to
# represent the gene -- chosen based on something like the mean absolute deviation
# (MAD) or the interquartile range (IQR). The other probe sets mapping to that 
# gene are then be disregarded and so information is lost. Alternatively,
# multiple probe sets can be combined into a new probe set using
# various summary statistics such as the mean or median. However,
# this approach also has caveats. The correct choice will depend
# on the specific data and analyses.

# Accepts a matrix expMat. The rownames of expMat are the "right" IDs.
# "mapping" is a list of vectors mapping "left" IDs to "right" IDs. The
# names of the list elements are the "left" IDs.
# "type" says what type of collapsing there will be. Will only one probe
# be selected for each gene or will multiple genes be summarized
# into a new probe? Options: "select" or "summary"
collapseExprMatrix = function( expMat, mapping, func=IQR, return.probes=FALSE,
                               type="select") {
  
  # make redundant "right" IDs that repeat according to how many map
  # to each "left" ID
  right = as.vector(unlist(mapping))
  right2 = right[ right %in% rownames(expMat) ]
  left = rep( names(mapping), times=sapply(mapping,length))
  left2 = left[ right %in% rownames(expMat) ]
  
  if(type=="select") {
    
    # group the scores from the chosen function by "left" IDs
    rowScore = as.matrix(apply(expMat,1,func))
    redundant.scores = rowScore[right2,]
    grouped.scores = split(redundant.scores,left2)
    
    # choose the "right" IDs that scored the highest
    foo = function(x) names(x)[which.max(x)]
    right.selected = sapply(grouped.scores,foo)
    
    # return the new expression matrix
    left.expMat = expMat[right.selected,]
    rownames(left.expMat) = names(right.selected)
    
    if(return.probes==TRUE) {
      return( list( expMat=left.expMat, collapseMap=right.selected ) )
    } else {
      return(left.expMat) 
    }
    
  } else if(type=="summary") {
    
    redundant.mat = expMat[right2,]
    dimnames(redundant.mat) = NULL
    ee = split( data.frame(redundant.mat),
                factor(left2) )
    left.expMat = t(sapply(ee,function(x) apply(x,2,func)))
    colnames(left.expMat) = colnames(expMat)
    return(left.expMat)
  }
}


# Accepts a vector or one-column matrix. See CollapseExpMat() for
# more details
collapseExprVector = function( expVector, mapping, func=median ) {
  
  # get the row-level score from the function input
  rowScore = expMat = as.matrix(expVector, ncol=1)

  # make redundant "right" IDs that repeat according to how many map
  # to each "left" ID
  right = as.vector(unlist(mapping))
  right2 = right[ right %in% rownames(expMat) ]
  left = rep( names(mapping), times=sapply(mapping,length))
  left2 = left[ right %in% rownames(expMat) ]
  
  # make a redundant matrix for easy collapsing
  redundant.scores = rowScore[right2,]
  grouped.scores = split(redundant.scores,left2)
  
  # collapse the probes according to the "func" function
  gene.ExprVector = sapply(grouped.scores,func)
  
  return(gene.ExprVector)  
}

