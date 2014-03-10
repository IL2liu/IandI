library(biomaRt)
library(limma)
library(Biostrings)
load("./data/annotation.RData/MsAnn.RData")
source("./C-Annotations/TransitiveMapping.R")

################################################################################
#                         Ensembl to MGI mapping                               #
################################################################################

## Jackson labs
raw = read.csv("./data/annotation.files/MGI-ID-to-Ensembl-MGI.rpt",sep="\t",header=F)
MGI.to.Ensembl.jackson = matrix(as.matrix(raw[!is.na(raw[,6]),c(1,6)]),ncol=2)

## Biomart
ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
filter = "mgi_id" # affy_mouse430_2 filter
features = as.character(unique(raw[,1]))
attribute = "ensembl_gene_id"
bm = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                      values = features, mart = ensembl ) )
MGI.to.Ensembl.biomart = matrix(bm[!is.na(bm[,2]),],ncol=2)

## Combine these databases
MGI.to.Ensembl.dup = rbind( MGI.to.Ensembl.jackson, MGI.to.Ensembl.biomart )
MGI.to.Ensembl.dup = trimWhiteSpace(MGI.to.Ensembl.dup)
duplicated = duplicated(as.data.frame(MGI.to.Ensembl.dup))
MGI.to.Ensembl = MGI.to.Ensembl.dup[ !duplicated, ]
MGI.to.Ensembl = MGI.to.Ensembl[ rowSums( MGI.to.Ensembl != "" ) == 2, ]
MGI.to.Ensembl[,1] = gsub("MGI:","",MGI.to.Ensembl[,1])

# save the annotation lists for MGI <-> Ensembl
MsAnn[["MGI"]][["Ensembl"]] = split(as.character(MGI.to.Ensembl[,2]),MGI.to.Ensembl[,1])
MsAnn[["Ensembl"]][["MGI"]] = inverseList(MsAnn[["MGI"]][["Ensembl"]])

# save the annotation lists for many-to-one mappings
# i.e. many Ensembl IDs to one MGI ID
summary(factor(sapply(MsAnn$Ensembl$MGI,length))) # MGI IDs/Ensembl ID "histogram"
summary(factor(sapply(MsAnn$MGI$Ensembl,length))) # Ensembl IDs/MGI ID "histogram"
# Will eventually need to map Ensembl probe sets to MGI probe sets.
# m2o = MsAnn$Ensembl$MGI[ sapply(MsAnn$Ensembl$MGI,length)==1 ] # only 
# MsAnn[["MGI"]][["Ensembl.many2one"]] = inverseList(m2o)
# MsAnn[["Ensembl.many2one"]][["MGI"]] = m2o

save(MsAnn,file="./data/annotation.RData/MsAnn.RData")