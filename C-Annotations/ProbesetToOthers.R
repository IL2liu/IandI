library(Matrix)
library(affy)
library(biomaRt)
load("./data/RData/esets.RData")
eset = esets[[1]]
library("mouse4302.db")
library(limma)


# Finding and downloading the DAVID files is not always that simple. Therefore,
# I uploaded every Affy ID to the web tool, performed the conversions, and
# downloaded the results.
write.table(featureNames(eset), 
            file="./all.affy.ids.txt",
            row.names=F, col.names=F, quote=F )

AffytoMGI.david = as.matrix(read.delim("./data/annotation.files/Affy-to-MGI-ID-DAVID.txt")[,1:2])
AffytoName.david = as.matrix(read.delim("./data/annotation.files/Affy-to-MGI-ID-DAVID.txt")[,c(1,4)])
AffytoSymbol.david = as.matrix(read.delim("./data/annotation.files/Affy-to-MGI-Symbol-DAVID.txt")[,1:2])
AffytoENTREZ.david = as.matrix(read.delim("./data/annotation.files/Affy-to-ENTREZ-ID-DAVID.txt")[,1:2])

# Biomart
# -------
# A package was written by Steffen Durinck and Wolgang Huber to easily access
# the biomaRt portal within R
ensembl = useMart("ensembl", dataset="mmusculus_gene_ensembl")
filter = "affy_mouse430_2" # affy_mouse430_2 filter

# Get a list of all the Affymetrix identifiers
features = featureNames(eset)

# Use the getBM function to query the database for the conversions
attribute = "mgi_symbol"
AffytoSymbol.biomart = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                                        values = features, mart = ensembl ) )
attribute = "mgi_description"
AffytoName.biomart = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                                      values = features, mart = ensembl ) )
attribute = "mgi_id"
AffytoMGI.biomart = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                                     values = features, mart = ensembl ) )
AffytoMGI.biomart[,2] = gsub("MGI:","",AffytoMGI.biomart[,2])
attribute = "entrezgene"
AffytoENTREZ.biomart = as.matrix( getBM(attributes = c(filter,attribute), filters = filter,
                                        values = features, mart = ensembl ) )



# Bioconductor
# ------------
# Bioconductor mappings
affy_ann = as.character( mouse4302SYMBOL )
AffytoSymbol.bioconductor = cbind( names(affy_ann), affy_ann  )
affy_ann = as.character( mouse4302GENENAME )
AffytoName.bioconductor = cbind( names(affy_ann), affy_ann )
affy_ann = as.character( mouse4302MGI )
AffytoMGI.bioconductor = cbind( names(affy_ann), affy_ann )
AffytoMGI.bioconductor[,2] = gsub("MGI:","",AffytoMGI.bioconductor[,2])
affy_ann = as.character( mouse4302ENTREZID )
AffytoENTREZ.bioconductor = cbind( names(affy_ann), affy_ann )


# Combine the three databases
# ---------------------------
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

colnames( AffytoMGI.david ) = c("From","To")
colnames( AffytoMGI.biomart ) = c("From","To")
colnames( AffytoMGI.bioconductor ) = c("From","To")
MGI = rbind( AffytoMGI.bioconductor, AffytoMGI.biomart, 
             AffytoMGI.david )

colnames( AffytoENTREZ.david ) = c("From","To")
colnames( AffytoENTREZ.biomart ) = c("From","To")
colnames( AffytoENTREZ.bioconductor ) = c("From","To")
ENTREZ = rbind( AffytoENTREZ.bioconductor, AffytoENTREZ.biomart, 
                AffytoENTREZ.david )

# remove any leading or trailing whitespace
Symbol = cbind( trimWhiteSpace(Symbol[,1]), trimWhiteSpace(Symbol[,2]) )
Name = cbind( trimWhiteSpace(Name[,1]), trimWhiteSpace(Name[,2]) )
MGI = cbind( trimWhiteSpace(MGI[,1]), trimWhiteSpace(MGI[,2]) )
ENTREZ = cbind( trimWhiteSpace(ENTREZ[,1]), trimWhiteSpace(ENTREZ[,2]) )

# remove duplicate rows
duplicated = duplicated(as.data.frame(Symbol))
Symbol = Symbol[ !duplicated, ]
duplicated = duplicated(as.data.frame(Name))
Name = Name[ !duplicated, ]
duplicated = duplicated(as.data.frame(MGI))
MGI = MGI[ !duplicated, ]
duplicated = duplicated(as.data.frame(ENTREZ))
ENTREZ = ENTREZ[ !duplicated, ]

# remove rows that have an empty cell
Name = Name[ rowSums( Name != "" ) == 2, ]
Symbol = Symbol[ rowSums( Symbol != "" ) == 2, ]
MGI = MGI[ rowSums( MGI != "" ) == 2, ]
ENTREZ = ENTREZ[ rowSums( ENTREZ != "" ) == 2, ]

# Put the annotations into list format
MsAnn = list()
MsAnn[["Affy"]][["Symbol"]] = split(as.character(Symbol[,2]),Symbol[,1])
MsAnn[["Affy"]][["Name"]] = split(as.character(Name[,2]),Name[,1])
MsAnn[["Affy"]][["MGI"]] = split(as.character(MGI[,2]),MGI[,1])
MsAnn[["Affy"]][["ENTREZ"]] = split(as.character(ENTREZ[,2]),ENTREZ[,1])
MsAnn[["Symbol"]][["Affy"]] = split(as.character(Symbol[,1]),Symbol[,2])
MsAnn[["Name"]][["Affy"]] = split(as.character(Name[,1]),Name[,2])
MsAnn[["MGI"]][["Affy"]] = split(as.character(MGI[,1]),MGI[,2])
MsAnn[["ENTREZ"]][["Affy"]] = split(as.character(ENTREZ[,1]),ENTREZ[,2])


# These annotations are stored in an R file
save(MsAnn,file = "./data/annotation.RData/MsAnn.RData")


