################################################################################
#                            Load data                                         #
################################################################################
library(affy)

data.directory = "./InVitro/data/CEL.files"
cel.files = grep("[.]CEL",list.files(data.directory,full.names=TRUE),value=TRUE)
abatch = ReadAffy(filenames=cel.files)

array.ids = gsub("[.]CEL","",sapply( strsplit(cel.files,"[/]"), function(x) x[length(x)] ))
phys.data = read.csv("./InVitro/data/CEL.files/fileDescriptions.csv",sep="\t")
pData(abatch) = phys.data[match(array.ids,phys.data$ID),c("Toxin","Time","Run")]
rownames(pData(abatch)) = array.ids
rownames(protocolData(abatch)@data) = rownames(phenoData(abatch)@data)

times = pData(abatch)$Time

save(abatch, file="./InVitro/data/RData/abatch.RData")