# Do preprocessing on each separate grouping of arrays
load("./data/RData/abatches.RData")
load("./data/annotation.RData/CustomPS.RData")
esets = list()
library(affy)
library(affyPLM)

################################################################################
#               Create CDF environment from annotation list                    #
################################################################################

#install.packages("/home/kevin/Downloads/mouse4302mmensgcdf_17.0.0.tar.gz")
library(mouse4302mmensgcdf)
ms = as.list(mouse4302mmensgcdf)
names(ms) = gsub("_at","",names(ms))
m4302 = list2env(ms)

################################################################################
#                           Preprocessing steps                                #
################################################################################


for( i in seq_along(abatches) ) {
  
  abatch = abatches[[i]]
  abatch@cdfName = "m4302"
  
  # background correction
  eset.bgc = bg.correct.rma(abatch)
  
  # probe normalization
  eset.norm = normalize.AffyBatch.loess(eset.bgc,type="pmonly",verbose=TRUE)

  # probe set summarization -- median polish
  eset.pss = rma(eset.norm, background=FALSE, normalize=FALSE)
  
  # probe set normalization
  eset.psn = normalize.ExpressionSet.loess(eset.pss)
  
  esets[[names(abatches)[i]]] = eset.psn
}

save(esets, file="./DiffProbeSets/data/RData/esets.RData")


