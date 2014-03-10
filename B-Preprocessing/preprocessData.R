# Do preprocessing on each separate grouping of arrays
load("./data/RData/abatches.RData")
esets = list()
library(affy)
library(affyPLM)

for( i in seq_along(abatches) ) {
  
  eset = abatches[[i]]
  
  # background correction
  eset.bgc = bg.correct.rma(eset)
  
  # probe normalization
  eset.norm = normalize.AffyBatch.loess(eset.bgc,type="pmonly",verbose=TRUE)

  # probe set summarization -- median polish
  eset.pss = rma(eset.norm, background=FALSE, normalize=FALSE)
  
  # probe set normalization
  eset.psn = normalize.ExpressionSet.loess(eset.pss)
  
  esets[[names(abatches)[i]]] = eset.psn
}

save(esets, file="./data/RData/esets.RData")


