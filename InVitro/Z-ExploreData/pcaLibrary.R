
################################################################################
#                               Functions                                      #
################################################################################

pca.fold = function( fold, pc.x=1, pc.y=2 ) {
  
  require("RColorBrewer")
  groups = factor(colnames(fold))
  centered = fold - rowMeans(fold)
  selected = rowSums( centered == 0 ) != ncol(fold)
  fold.pca = fold[selected,]
  
  pca = prcomp( t(fold.pca), center=TRUE, scale=TRUE )
  pca.info = summary(pca)
  
  num.groups = nlevels(groups)
  if( num.groups < 9 ) {
    point.color.palette = brewer.pal(num.groups,"Dark2")
  } else {
    point.color.palette = rep(brewer.pal(8,"Dark2"),100)[1:num.groups]
  }
  point.colors = point.color.palette[ as.numeric(groups) ]
  point.shapes = as.numeric(groups)
  
  # the plot
  layout(matrix(1:2, nc = 2), c(3,1))
  plot(pca$x[,c(pc.x,pc.y)],
       main="Principal component analysis", col=point.colors,
       xlab=paste("PC1:", round(pca.info$imp[2,pc.x]*100, digits=2), "%"),
       ylab=paste("PC2:", round(pca.info$imp[2,pc.y]*100, digits=2), "%"),
       pch=point.shapes, lwd=5, cex=3)
  xmax = max(pca$x[,pc.x])
  ymax = max(pca$x[,pc.y])
  xrange = abs(xmax - min(pca$x[,pc.x]))
  par(xpd=NA)
  legend(xmax+0.1*xrange,ymax, legend=levels(groups), 
         col=brewer.pal(num.groups,"Dark2"),
         pch = 1:num.groups,
         bty="n", cex=2, pt.lwd=12)
  
}


# A very specific function for performing PCA and plotting results
pca.plot = function( eset, pc.x=1, pc.y=2 ) {
  
  require("RColorBrewer")
  groups = interaction( pData(eset)$Toxin, pData(eset)$Time )
  
  # I want to scale the data, but genes that are the exact same number can not
  # be scaled to unit variance. Hence, I will remove these variables
  centered = exprs(eset) - rowMeans(exprs(eset))
  selected = rowSums( centered == 0 ) != ncol(eset)
  eset.pca = eset[selected,]
  
  # Now, actually run the PCA function
  pca = prcomp( t(exprs(eset.pca)), center=TRUE, scale=TRUE )
  pca.info = summary(pca)
  
  # Set up the options for plotting the points
  num.groups = nlevels(groups)
  if( num.groups < 9 ) {
    point.color.palette = brewer.pal(num.groups,"Dark2")
  } else {
    point.color.palette = rep(brewer.pal(8,"Dark2"),100)[1:num.groups]
  }
  point.colors = point.color.palette[ as.numeric(groups) ]
  point.shapes = as.numeric(groups)
  
  # if i want to color by the run
  point.colors = point.color.palette[ pData(eset)$Run  ]
  
  # the plot
  layout(matrix(1:2, nc = 2), c(3,1))
  plot(pca$x[,c(pc.x,pc.y)],
       main="Principal component analysis", col=point.colors,
       xlab=paste("PC1:", round(pca.info$imp[2,pc.x]*100, digits=2), "%"),
       ylab=paste("PC2:", round(pca.info$imp[2,pc.y]*100, digits=2), "%"),
       pch=point.shapes, lwd=5, cex=3)
  xmax = max(pca$x[,pc.x])
  ymax = max(pca$x[,pc.y])
  xrange = abs(xmax - min(pca$x[,pc.x]))
  par(xpd=NA)
  legend(xmax+0.1*xrange,ymax, legend=levels(groups), 
         col=brewer.pal(num.groups,"Dark2"),
         pch = 1:num.groups,
         bty="n", cex=2, pt.lwd=6)
}






plot.fit = function( eset, fit ) {
  newset = exprs(eset)
  for( i in 1:nrow(exprs(eset))) {
    resid = exprs(eset)[i,] - fit$design %*% fit$coefficients[i,]
    newset[i,] =  fit$design[,1:9] %*% fit$coefficients[i,1:9] + resid
  }
  eset2 = eset
  exprs(eset2) = newset
  pca.plot(eset2, 1, 2)  
}
