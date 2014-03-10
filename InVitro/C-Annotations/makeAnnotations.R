library(Matrix)
library(affy)
library(biomaRt)
load("./InVitro/data/RData/bppsl.eset.RData")
eset = bppsl.eset
library("hgu133plus2.db")
library(limma)

################################################################################
#                         Gene IDs                                             #
################################################################################

# Finding and downloading the DAVID files is not always that simple. Therefore,
# I uploaded every Affy ID to the web tool, performed the conversions, and
# downloaded the results.
write.table(featureNames(eset), 
            file="./InVitro/data/annotation.files/all.affy.ids.txt",
            row.names=F, col.names=F, quote=F )

AffytoSymbol.david = as.matrix(read.delim("./InVitro/data/annotation.files/Affy-to-HGNC-DAVID.txt")[,1:2])
AffytoName.david = as.matrix(read.delim("./InVitro/data/annotation.files/Affy-to-HGNC-DAVID.txt")[,c(1,4)])
AffytoENTREZ.david = as.matrix(read.delim("./InVitro/data/annotation.files/Affy-to-ENTREZ-human-DAVID.txt")[,1:2])


################ end here so far
# Biomart
# A package was written by Steffen Durinck and Wolgang Huber to easily access
# the biomaRt portal within R
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
filter = "affy_hg_u133_plus_2" # affy_mouse430_2 filter

# Get a list of all the Affymetrix identifiers
features = featureNames(eset)

# Use the getBM function to query the database for the conversions
attribute = "hgnc_symbol"
AffytoSymbol.biomart = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                                        values = features, mart = ensembl ) )
attribute = "description"
AffytoName.biomart = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                                      values = features, mart = ensembl ) )
attribute = "entrezgene"
AffytoENTREZ.biomart = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                                        values = features, mart = ensembl ) )

# Bioconductor
# Bioconductor mappings
affy_ann = as.character( hgu133plus2SYMBOL )
AffytoSymbol.bioconductor = cbind( names(affy_ann), affy_ann  )
affy_ann = as.character( hgu133plus2GENENAME )
AffytoName.bioconductor = cbind( names(affy_ann), affy_ann )
affy_ann = as.character( hgu133plus2ENTREZID )
AffytoENTREZ.bioconductor = cbind( names(affy_ann), affy_ann )


# Combine the three databases
colnames( AffytoSymbol.david ) = c("From","To")
colnames( AffytoSymbol.biomart ) = c("From","To")
colnames( AffytoSymbol.bioconductor ) = c("From","To")
Symbol = rbind( AffytoSymbol.bioconductor, AffytoSymbol.biomart, 
                AffytoSymbol.david )

colnames( AffytoName.david ) = c("From","To")
colnames( AffytoName.biomart ) = c("From","To")
colnames( AffytoName.bioconductor ) = c("From","To")
Name = rbind( AffytoName.bioconductor, AffytoName.biomart, 
              AffytoName.david )

colnames( AffytoENTREZ.david ) = c("From","To")
colnames( AffytoENTREZ.biomart ) = c("From","To")
colnames( AffytoENTREZ.bioconductor ) = c("From","To")
ENTREZ = rbind( AffytoENTREZ.bioconductor, AffytoENTREZ.biomart, 
                AffytoENTREZ.david )

# remove any leading or trailing whitespace
Symbol = cbind( trimWhiteSpace(Symbol[,1]), trimWhiteSpace(Symbol[,2]) )
Name = cbind( trimWhiteSpace(Name[,1]), trimWhiteSpace(Name[,2]) )
ENTREZ = cbind( trimWhiteSpace(ENTREZ[,1]), trimWhiteSpace(ENTREZ[,2]) )

# Move all symbols to uppercase
Symbol[,2] = toupper(Symbol[,2])

# remove duplicate rows
duplicated = duplicated(as.data.frame(Symbol))
Symbol = Symbol[ !duplicated, ]
duplicated = duplicated(as.data.frame(Name))
Name = Name[ !duplicated, ]
duplicated = duplicated(as.data.frame(ENTREZ))
ENTREZ = ENTREZ[ !duplicated, ]

# remove rows that have an empty cell
Name = Name[ rowSums( Name != "" ) == 2, ]
Symbol = Symbol[ rowSums( Symbol != "" ) == 2, ]
ENTREZ = ENTREZ[ rowSums( ENTREZ != "" ) == 2, ]

# Put the annotations into list format
HsAnn = list()
HsAnn[["Affy"]][["Symbol"]] = split(as.character(Symbol[,2]),Symbol[,1])
HsAnn[["Affy"]][["Name"]] = split(as.character(Name[,2]),Name[,1])
HsAnn[["Affy"]][["ENTREZ"]] = split(as.character(ENTREZ[,2]),ENTREZ[,1])
HsAnn[["Symbol"]][["Affy"]] = split(as.character(Symbol[,1]),Symbol[,2])
HsAnn[["Name"]][["Affy"]] = split(as.character(Name[,1]),Name[,2])
HsAnn[["ENTREZ"]][["Affy"]] = split(as.character(ENTREZ[,1]),ENTREZ[,2])

# These annotations are stored in an R file
save(HsAnn,file = "./InVitro/data/annotation.RData/HumanAnnotations.RData")

################################################################################
#                       Gene ontologies                                        #
################################################################################

# The GO Annotations are from the July 11, 2012 update of the Gene Ontology
aa = read.delim("./InVitro/data/annotation.files/gene_association.goa_human",header=F)
bb = aa[,c(3,4,5,7,9)]
colnames(bb) = c("HGNC","NOT?","GOid","Evidence","Ontology")

# Remove NOT rows
cc = bb[!(bb[,2] %in% c("NOT","NOT|contributes_to")),]

# Remove rows with ND evidence code
cc = cc[ cc$Evidence != "ND", ]

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
HsAnnGO = list()
HsAnnGO[["Symbol"]][["BP"]] = split(as.character(BP.nondup[,2]),BP.nondup[,1])
HsAnnGO[["BP"]][["Symbol"]] = split(as.character(BP.nondup[,1]),BP.nondup[,2])
HsAnnGO[["Symbol"]][["MF"]] = split(as.character(MF.nondup[,2]),MF.nondup[,1])
HsAnnGO[["MF"]][["Symbol"]] = split(as.character(MF.nondup[,1]),MF.nondup[,2])
HsAnnGO[["Symbol"]][["CC"]] = split(as.character(CC.nondup[,2]),CC.nondup[,1])
HsAnnGO[["CC"]][["Symbol"]] = split(as.character(CC.nondup[,1]),CC.nondup[,2])

# Also make lists for non-IEA annotations
ee = dd[dd$Evidence!="IEA",]
BP = as.matrix(ee[ee$Ontology == "P",1:3])
MF = as.matrix(ee[ee$Ontology == "F",1:3])
CC = as.matrix(ee[ee$Ontology == "C",1:3])
BP.nondup = BP[!duplicated(BP[,1:2]),]
MF.nondup = MF[!duplicated(MF[,1:2]),]
CC.nondup = CC[!duplicated(CC[,1:2]),]
HsAnnGO[["Symbol"]][["BP.noIEA"]] = split(as.character(BP.nondup[,2]),BP.nondup[,1])
HsAnnGO[["BP.noIEA"]][["Symbol"]] = split(as.character(BP.nondup[,1]),BP.nondup[,2])
HsAnnGO[["Symbol"]][["MF.noIEA"]] = split(as.character(MF.nondup[,2]),MF.nondup[,1])
HsAnnGO[["MF.noIEA"]][["Symbol"]] = split(as.character(MF.nondup[,1]),MF.nondup[,2])
HsAnnGO[["Symbol"]][["CC.noIEA"]] = split(as.character(CC.nondup[,2]),CC.nondup[,1])
HsAnnGO[["CC.noIEA"]][["Symbol"]] = split(as.character(CC.nondup[,1]),CC.nondup[,2])

save(HsAnnGO,file = "./InVitro/data/annotation.RData/HumanGOAnnotations.RData")
remove(list=ls())

