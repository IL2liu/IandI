################################################################################
#                     MGI to symbol mapping                                    #
################################################################################

load("./data/annotation.RData/MsAnn.RData")
library(mouse4302.db)
library(biomaRt)
library(limma)

## Jackson labs
raw = read.csv("./data/annotation.files/MGI-to-Gene-Symbol-MGI.rpt",sep="\t",header=T)
MGI2Symbol.jackson = as.matrix(raw[!is.na(raw[,7]),c(1,7)])

## Biomart
ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
filter = "mgi_id" # affy_mouse430_2 filter
features = as.character(unique(raw[,1]))
attribute = "mgi_symbol"
bm = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                      values = features, mart = ensembl ) )
MGI2Symbol.biomart = bm[!is.na(bm[,2]),]

## Combine the three datbases
colnames( MGI2Symbol.jackson ) = c("From","To")
colnames( MGI2Symbol.biomart ) = c("From","To")
MGI2EG.all = rbind( MGI2Symbol.jackson, MGI2Symbol.biomart )
MGI2EG.all = trimWhiteSpace(MGI2EG.all)
duplicated = duplicated(as.data.frame(MGI2EG.all))
MGI2EG = MGI2EG.all[ !duplicated, ]
MGI2EG = MGI2EG[ rowSums( MGI2EG != "" ) == 2, ]
MGI2EG[,1] = gsub("MGI:","",MGI2EG[,1])
MsAnn[["MGI"]][["Symbol"]] = split(as.character(MGI2EG[,2]),MGI2EG[,1])
MsAnn[["Symbol"]][["MGI"]] = split(as.character(MGI2EG[,1]),MGI2EG[,2])

save(MsAnn,file="./data/annotation.RData/MsAnn.RData")





