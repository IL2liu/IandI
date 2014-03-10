################################################################################
#                 Change the format of gene sets for CAMERA                    #
################################################################################
indexGeneSets = function(genesets, collapseMap, minSize=10, maxSize=700 ) {
  
  # CAMERA needs the gene sets to be made by the indeces of the test statistic
  # vector. For example, if gene A is the fifth gene in the in the vector
  # input to CAMERA, then all gene A's should be replaced with "5" in the 
  # geneset annotation list
  gene.names = names(collapseMap)
  indeces = match(unlist(genesets),gene.names)
  redundant.gene.set.names = rep( names(genesets), times=sapply(genesets,length))
  mapping = na.omit( cbind(redundant.gene.set.names,indeces) )
  genesets.indeces = split(as.numeric(mapping[,2]),mapping[,1])
  
  # Next, it is possible that one probe set represents many
  # genes in a geneset. This may have happened while
  # an expression matrix was collapsed from probe sets to gene IDs.
  # This is corrected in the code below.
  foo = function(x) {
    probesets = as.character(collapseMap[x])
    # duplicated function keeps the first probe set
    return(x[ !duplicated(probesets) ]) 
  }
  genesets.indeces = lapply(genesets.indeces,foo)

  # Set the limit on the smallest and largest allowable gene sets
  geneSets.cut = genesets.indeces[ sapply(genesets.indeces,length)>=minSize ]
  geneSets.cut = geneSets.cut[ sapply(geneSets.cut,length)<=maxSize ]
  return(geneSets.cut)
}

