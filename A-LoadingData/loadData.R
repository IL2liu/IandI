library(affy)

# load in the expression data 
data.directory = "./data/CEL.files"
cel.files = grep("[.]CEL",list.files(data.directory,full.names=TRUE),value=TRUE)
abatch = ReadAffy(filenames=cel.files)

# and other data associated with each sample
array.ids = gsub("[.]CEL","",gsub(".*_","",cel.files))
phys.data = read.csv("./PhysiologyMeasurements.csv",sep="\t")
pData(abatch) = phys.data[match(array.ids,phys.data$ID),c("Toxin","Time")]
rownames(pData(abatch)) = array.ids
rownames(protocolData(abatch)@data) = rownames(phenoData(abatch)@data)

# For when comparisons are being made to identify differentially expressed genes,
# each time point is considered separately. 
times = pData(abatch)$Time
abatches = list()
abatches[["2hr"]] = abatch[,times==2]
abatches[["6hr"]] = abatch[,times==6]
abatches[["16hr"]] = abatch[,times==16]
abatches[["all"]] = abatch

save(abatches,file="./data/RData/abatches.RData")

