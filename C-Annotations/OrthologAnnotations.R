library(limma)
source("./C-Annotations/TransitiveMapping.R")
load("./data/annotation.RData/MsAnn.RData")

################################################################################
#              Mouse ENTREZ to Human ENTREZ from NCBI Homologene               #
################################################################################

# Human Taxonomy ID: 9606
# Mus Musculus ID: 10090
aa = trimWhiteSpace(as.matrix(read.delim(
  file="./data/annotation.files/Mouse-ENTREZ-to-Human-ENTREZ-NCBI-homologene.txt",header=F)))
colnames(aa) = c("HomologID","TaxonomyID","EntrezID","Symbol","Protein gi","ProteinAccession")
human = aa[aa[,2]=="9606",c(1,3)] # HomologID -> Human Entrez ID
mouse = aa[aa[,2]=="10090",c(1,3)] # HomologID -> Mouse Entrez ID

# Now do a transitive mapping: Mouse Entrez ID -> HomologID -> Human Entrez ID
mouse2human.ENTREZ.homologene = TransitiveMapping(mouse,human)



################################################################################
#               Mouse ENTREZ to Human ENTREZ from MGI file                     #
################################################################################

aa = trimWhiteSpace(as.matrix(read.delim(
          file="./data/annotation.files/Mouse-ENTREZ-to-Human-ENTREZ-MGI.txt",header=T)))
mouse2human.ENTREZ.jackson = aa[,c(7,2)]



################################################################################
#                 Combine NCBI and MGI ortholog annotations                    #
################################################################################
M2H = rbind(mouse2human.ENTREZ.homologene,
            mouse2human.ENTREZ.jackson)
duplicated = duplicated(as.data.frame(M2H))
M2H = M2H[ !duplicated, ]
M2H = M2H[ rowSums(is.na(M2H))==0, ]
Mm.to.Hs.ENTREZ = M2H[ rowSums( M2H != "" ) == 2, ] #remove empty rows

dd = Mm.to.Hs.ENTREZ
MsAnn[["HumanENTREZ"]][["ENTREZ"]] = split(dd[,1],dd[,2])
MsAnn[["ENTREZ"]][["HumanENTREZ"]] = split(dd[,2],dd[,1])


################################################################################
#                    Affy probe set to Human ENTREZ                            #
################################################################################

# Affy probe set -> ENTREZ -> Human ENTREZ
dd = TransitiveMapping( MsAnn$ENTREZ$Affy, MsAnn$ENTREZ$HumanENTREZ  )
MsAnn[["Affy"]][["HumanENTREZ"]] = dd
MsAnn[["HumanENTREZ"]][["Affy"]] = inverseList(dd)


save(MsAnn,file = "./data/annotation.RData/MsAnn.RData")


################################################################################
#                       ENSEMBL Human-Mouse Orthologs                          #
################################################################################

# orthology mappings
OrthoAnn = list()
aa = read.csv("./data/annotation.files/Human-ENSG-to-Mouse-ENSMUSG.txt",
              stringsAsFactors=FALSE)
bb = aa[,c(1,3)]
duplicated = duplicated(as.data.frame(bb))
cc = bb[ !duplicated, ]
dd = cc

# save the mappings
OrthoAnn$hs$mm = split(as.character(dd[,2]),dd[,1])
OrthoAnn$mm$hs = split(as.character(dd[,1]),dd[,2])

# Mouse to Human. One-to-one
mm.to.one = names(OrthoAnn$mm$hs)[sapply(OrthoAnn$mm$hs, length)==1]
one.to.hs = OrthoAnn$hs$mm[ unlist(OrthoAnn$mm$hs[mm.to.one]) ]
hs.to.mm.one2one = one.to.hs[ sapply(one.to.hs,length)==1 ]
one2one = as.list( names(hs.to.mm.one2one) )
names(one2one) = unlist(hs.to.mm.one2one)

# make maps excluding the one to one mappings
hs.mm = OrthoAnn$hs$mm[ !(names(OrthoAnn$hs$mm) %in% unlist(one2one)) ]
mm.hs = OrthoAnn$mm$hs[ !(names(OrthoAnn$mm$hs) %in% names(one2one)) ]

# next find one to many (mouse to human)
tomany = names(mm.hs[ sapply( mm.hs, length ) != 1 ])
rtomany = names(hs.mm[ sapply( hs.mm, length ) != 1 ])
one2many = mm.hs[ names(which(sapply( mm.hs, function(x) !all(x %in% rtomany) ))) ]

# next find many to one
many2one = hs.mm[ names(which(sapply(hs.mm, function(x) !all(x %in% tomany) ))) ]

# naturally, the rest should be many to many relationships
others = c(names(one2one),names(one2many),unique(unlist(many2one)))
m2mps = setdiff( names(mm.hs), others )
many2many = mm.hs[ m2mps ]

save(OrthoAnn, file="./data/annotation.RData/OrthoAnn.RData")

