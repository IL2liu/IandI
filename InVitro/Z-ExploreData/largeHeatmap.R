# show the top genes from different analyses
heatmap.ibmt.cor = function( res ) {
  require(gplots)
  lfc = res$coefficients
  
  distfun = function(x) as.dist( 1 - cor(x) )
  hclustfun = function(x) hclust(x,method="ward")
  
  # correlation plot
  heatmap.2( cor(lfc,method="pearson"),
             col=colorpanel(100,"blue","white","maroon"),
             symbreaks=TRUE,
             trace="none",
             distfun=distfun,
             hclustfun=hclustfun,
             cellnote=formatC(cor(lfc),format="f",digits=2,width=1),
             notecol="black",
             density.info="none", dendrogram="none", Rowv=FALSE, Colv=FALSE)
}


# show the top genes from different analyses
heatmap.ibmt.genes = function( res ) {
  require(gplots)
  pp = res$IBMT.p
  ap = apply(pp,2,function(x) p.adjust(x,"BH"))
  de.list = apply( ap < 0.01, 2, function(x) names(which(x)))
  lfc = res$coefficients
  
  allgenes = unique(unlist(de.list))
  exp.mat = lfc[allgenes,]
  
  distfun = function(x) as.dist( 1 - cor(t(x)) )
  hclustfun = function(x) hclust(x,method="ward")
  
  # correlation plot
  heatmap.2( exp.mat,
             col=colorpanel(100,"blue","white","red"),
             symbreaks=TRUE,
             trace="none",
             distfun=distfun,
             hclustfun=hclustfun,
             density.info="none", Colv=FALSE, dendrogram="none",
             labRow=NA)
}



