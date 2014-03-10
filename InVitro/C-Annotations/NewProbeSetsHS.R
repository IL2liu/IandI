
library(biomaRt)
library(limma)
library(Biostrings)
load("./InVitro/data/annotation.RData/HumanAnnotations.RData")
source("./C-Annotations/TransitiveMapping.R")
library(hgu133plus2cdf)
library(hgu133plus2probe)
library(affy)

################################################################################
#    Ensembl transcript to Ensembl Gene (from Ensembl's database file)         #
################################################################################

# Ensembl transcript versus gene IDs
Ensembl.transcript.to.gene = read.table("./InVitro/data/annotation.files/Ensembl-transcript-IDs-to-Ensembl-gene-IDs-Ensembl.txt")
sum(duplicated(Ensembl.transcript.to.gene[,1])) # shows that every transcript
# ID is unique in the transcript to gene mapping. Hence, the transcript
# to gene mapping is many-to-one. In other words, each transcript ID
# maps to only one gene. Thus, the probes for all of a gene's transcripts 
# can then represent the expression of that gene.
#  sed -rn 's/.*ENST([0-9]*).*ENSG([0-9]*).*/ENST\1\tENSG\2/p' 
#  ./Homo_sapiens.GRCh37.71.cdna.all.fa >Ensembl-transcript-IDs-to-Ensembl-gene-IDs-Ensembl.txt
# This SED command on the fasta file below generates the text mapping

################################################################################
#                   Getting the sequences for each probe                       #
################################################################################


# get the probe sequences from bioconductor
probes = as.data.frame(hgu133plus2probe)

# rename the probe sequences by their probe ID (not x,y position)
positions = as.matrix(probes[,c("x","y")])
ids = apply(positions,1,function(x) xy2indices(x[1],x[2],cdf="hgu133plus2cdf"))
probe.seqs = probes[,"sequence"]
names(probe.seqs) = ids

# Prepare fasta file for BLAST search
writeXStringSet(DNAStringSet(probe.seqs),"./InVitro/data/annotation.files/ProbeSeqs.fasta")

################################################################################
#                 BLAST probe sequences to Ensembl transcripts                 #
################################################################################

# The blast executables and databases are far too large to include
# in a download we can provide. However, here are the instructions
# for repeating what we did.
# ncbi-blast-2.2.28+-x64-linux.tar.gz from the NCBI blast website was unzipped to 
# folder /anypath/
# Mus_musculus_.GRCm38.70.cdna.all.fa.gz from the ENSEMBL website was downloaded 
# to /anypath/db/ , extracted, and used to create a blast database.
# ./bin/makeblastdb -dbtype 'nucl' -in ./db/Homo_sapiens.GRCh37.71.cdna.all.fa 
# -out ensembl -input_type fasta -title ensembl -hash_index -taxid 9606

# /anypath/bin/makblastdb -dbtype 'nucl' -in /anypath/db/Mus_musculus_.GRCm38.70.cnda.all.fa
#                         -out ensembl -input_type fasta -title ensembl -hash_index -taxid 10090
# The ProbeSeqs.fasta file created above was copied to /anypath/db/ProbeSeqs.fasta

# ./blastn -evalue 10000 -word_size 25 -perc_identity 100 -db ensembl 
#          -query ../db/ProbeSeqs.fasta -out ../db/mappings.txt -outfmt '6 qseqid sseqid'
# Then copy mappings.txt to the ./data/annotation.files/Probe-to-Ensembl-transcript-ID.txt


################################################################################
#           Parsing results from BLAST of probes to Ensembl transcripts        #
################################################################################


# convert the Ensembl transcripts from blast to Ensembl genes. This is
# simple since the mapping is many-to-one
probe.to.ensembl.transcript = 
  read.delim("./InVitro/data/annotation.files/AffyProbe-to-Ensembl-transcript-ID.txt",header=FALSE)
blast.transcripts = as.character(probe.to.ensembl.transcript[,2])
blast.genes = as.character( Ensembl.transcript.to.gene[ match(blast.transcripts, 
                                                              Ensembl.transcript.to.gene[,1] ), 2] )
probe.to.ensembl.gene = probe.to.ensembl.transcript
probe.to.ensembl.gene[,2] = blast.genes

# We now have a potentially many-to-many mapping from Affy probe to Ensembl gene.
# We reduce this to a many-to-many mapping by removing probes
probe.to.gene = probe.to.ensembl.gene[ !duplicated(probe.to.ensembl.gene), ]
probe.counts = table(probe.to.gene[,1])
many2one.probes = names(probe.counts)[probe.counts==1]
p2g.many2one  = probe.to.gene[ probe.to.gene[,1]  %in% many2one.probes, ]
probe.to.gene = na.omit(p2g.many2one)


################################################################################
#           Defining custom MGI or Ensembl gene base probe sets                #
################################################################################

CustomPS = list()

# Ensembl probe sets
CustomPS[["Ensembl"]][["probe"]] = split(probe.to.gene[,1],probe.to.gene[,2])
CustomPS[["probe"]][["Ensembl"]] = inverseList(CustomPS[["Ensembl"]][["probe"]])

save(CustomPS,file="./InVitro/data/annotation.RData/CustomPS.RData")


################################################################################
#                  Ensembl Gene IDs to ENTREZ Gene IDs                         #
################################################################################

ens = unique(as.character(Ensembl.transcript.to.gene[,2]))
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
attribute = "entrezgene"
filter = "ensembl_gene_id"

# Use the getBM function to query the database for the conversions
ens2ez.bm = getBM(filters=filter, attributes=c(filter,attribute), values = ens, mart = mart )
entrez.mart = ens2ez.bm

# Also use DAVID
#ens2 = names(CustomPS$Ensembl$probe)
#write.table(ens2,row.names=FALSE,col.names=FALSE,quote=FALSE,file="/home/kevin/test.txt")
#aa = read.delim("./InVitro/data/annotation.files/Ensembl-to-Entrez-human-DAVID.txt")
#entrez.david = aa[,1:2]

# Also use Bioconductor
require(org.Hs.eg.db)
aa = as.list(org.Hs.egENSEMBL)
bb = aa[ sapply(aa,function(x) !any(is.na(x)) ) ]
entrez.bc = cbind( unlist(bb), rep(names(bb),times=sapply(bb,length)) )

# Combine the three databases
colnames( entrez.mart ) = c("From","To")
#colnames( entrez.david ) = c("From","To")
colnames( entrez.bc ) = c("From","To")
entrez = rbind( entrez.mart, entrez.bc )

# Combining the databases
Entrez = cbind( trimWhiteSpace(entrez[,1]), trimWhiteSpace(entrez[,2]) )
duplicated = duplicated(as.data.frame(Entrez))
Entrez = Entrez[ !duplicated, ]
Entrez = Entrez[ !apply(is.na(Entrez), 1, any), ]
Entrez = Entrez[ rowSums( Entrez != "" ) == 2, ]
dimnames(Entrez) = NULL

r1 = names(CustomPS$Ensembl$probe)
r2 = Entrez[,1]
length(setdiff(r1,r2))
length(setdiff(r2,r1))
length(intersect(r1,r2))

# Now figure out what to do with Ensembl IDs that map to many Entrez IDs
aa = split(as.character(Ent[,2]),Ent[,1])
table(sapply(aa,length))
aa["ENSG00000161992"]
aa[sapply(aa, length)>1]

# Now save the results
HsAnn[["Ensembl"]][["Entrez"]] = split(as.character(Entrez[,2]),Entrez[,1])
HsAnn[["Entrez"]][["Ensembl"]] = split(as.character(Entrez[,2]),Entrez[,1])

################################################################################
#             Ensembl Gene IDs to Official Gene Symbols                        #
################################################################################

ens = unique(as.character(Ensembl.transcript.to.gene[,2]))
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
attribute = "hgnc_symbol"
filter = "ensembl_gene_id"

# Use the getBM function to query the database for the conversions
hgnc.bm = getBM(filters=filter, attributes=c(filter,attribute), values = ens, mart = mart )

HGNC = cbind( trimWhiteSpace(hgnc.bm[,1]), trimWhiteSpace(hgnc.bm[,2]) )
duplicated = duplicated(as.data.frame(HGNC))
HGNC = HGNC[ !duplicated, ]
HGNC = HGNC[ !apply(is.na(HGNC), 1, any), ]
HGNC = HGNC[ rowSums( HGNC != "" ) == 2, ]
dimnames(HGNC) = NULL

HsAnn[["Ensembl"]][["Symbol"]] = split(as.character(HGNC[,2]),HGNC[,1])
HsAnn[["Symbol"]][["Ensembl"]] = split(as.character(HGNC[,1]),HGNC[,2])

save(HsAnn,file="./InVitro/data/annotation.RData/HumanAnnotations.RData")












