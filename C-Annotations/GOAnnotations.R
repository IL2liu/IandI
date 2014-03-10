library(limma)

# The GO Annotations are from the July 11, 2012 update of the Gene Ontology
aa = read.delim("./data/annotation.files/GO-to-MGI-ID.mgi",header=F)
bb = aa[,c(2,4,5,7,9)]
colnames(bb) = c("MGIid","NOT?","GOid","Evidence","Ontology")

# Remove NOT rows
cc = bb[!(bb[,2] %in% c("NOT","NOT|contributes_to")),]

# Remove rows with ND evidence code
cc = cc[ cc$Evidence != "ND", ]

# Remove the characters "MGI:"
cc[,1] = trimWhiteSpace(gsub("MGI:","",cc[,1]))

# Remove the NOT row
dd = cc[,-2]

# separate the 3 ontologies
BP = as.matrix(dd[dd$Ontology == "P",1:3])
MF = as.matrix(dd[dd$Ontology == "F",1:3])
CC = as.matrix(dd[dd$Ontology == "C",1:3])

# remove redundant annotations made because an annotation has multiple
# evidence codes
BP.nondup = BP[!duplicated(BP[,1:2]),]
MF.nondup = MF[!duplicated(MF[,1:2]),]
CC.nondup = CC[!duplicated(CC[,1:2]),]

# Put the annotations into list format
MsAnnGO = list()
MsAnnGO[["MGI"]][["BP"]] = split(as.character(BP.nondup[,2]),BP.nondup[,1])
MsAnnGO[["BP"]][["MGI"]] = split(as.character(BP.nondup[,1]),BP.nondup[,2])
MsAnnGO[["MGI"]][["MF"]] = split(as.character(MF.nondup[,2]),MF.nondup[,1])
MsAnnGO[["MF"]][["MGI"]] = split(as.character(MF.nondup[,1]),MF.nondup[,2])
MsAnnGO[["MGI"]][["CC"]] = split(as.character(CC.nondup[,2]),CC.nondup[,1])
MsAnnGO[["CC"]][["MGI"]] = split(as.character(CC.nondup[,1]),CC.nondup[,2])

# Also make lists for non-IEA annotations
ee = dd[dd$Evidence!="IEA",]
BP = as.matrix(ee[ee$Ontology == "P",1:3])
MF = as.matrix(ee[ee$Ontology == "F",1:3])
CC = as.matrix(ee[ee$Ontology == "C",1:3])
BP.nondup = BP[!duplicated(BP[,1:2]),]
MF.nondup = MF[!duplicated(MF[,1:2]),]
CC.nondup = CC[!duplicated(CC[,1:2]),]
MsAnnGO[["MGI"]][["BP.noIEA"]] = split(as.character(BP.nondup[,2]),BP.nondup[,1])
MsAnnGO[["BP.noIEA"]][["MGI"]] = split(as.character(BP.nondup[,1]),BP.nondup[,2])
MsAnnGO[["MGI"]][["MF.noIEA"]] = split(as.character(MF.nondup[,2]),MF.nondup[,1])
MsAnnGO[["MF.noIEA"]][["MGI"]] = split(as.character(MF.nondup[,1]),MF.nondup[,2])
MsAnnGO[["MGI"]][["CC.noIEA"]] = split(as.character(CC.nondup[,2]),CC.nondup[,1])
MsAnnGO[["CC.noIEA"]][["MGI"]] = split(as.character(CC.nondup[,1]),CC.nondup[,2])


save(MsAnnGO,file = "./data/annotation.RData/MsAnnGO.RData")
remove(list=ls())
