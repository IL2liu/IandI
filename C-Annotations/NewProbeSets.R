
library(biomaRt)
library(limma)
library(Biostrings)
load("./data/annotation.RData/MsAnn.RData")
source("./C-Annotations/TransitiveMapping.R")
library(mouse4302cdf)
library(mouse4302probe)
library(affy)

################################################################################
#    Ensembl transcript to Ensembl Gene (from Ensembl's database file)         #
################################################################################

# Ensembl transcript versus gene IDs
Ensembl.transcript.to.gene = read.table("./data/annotation.files/Ensembl-transcript-IDs-to-Ensembl-gene-IDs-Ensembl.txt")
sum(duplicated(Ensembl.transcript.to.gene[,1])) # shows that every transcript
# ID is unique in the transcript to gene mapping. Hence, the transcript
# to gene mapping is many-to-one. In other words, each transcript ID
# maps to only one gene. Thus, the probes for all of a gene's transcripts 
# can then represent the expression of that gene.


################################################################################
#                   Getting the sequences for each probe                       #
################################################################################


# get the probe sequences from bioconductor
probes = as.data.frame(mouse4302probe)

# rename the probe sequences by their probe ID (not x,y position)
positions = as.matrix(probes[,c("x","y")])
ids = apply(positions,1,function(x) xy2indices(x[1],x[2],cdf="mouse4302cdf"))
probe.seqs = probes[,"sequence"]
names(probe.seqs) = ids

# Prepare fasta file for BLAST search
writeXStringSet(DNAStringSet(probe.seqs),"./data/annotation.files/ProbeSeqs.fasta")


################################################################################
#                 BLAST probe sequences to Ensembl transcripts                 #
################################################################################

# The blast executables and databases are far too large to include
# in a download we can provide. However, here are the instructions
# for repeating what we did.
# ncbi-blast-2.2.27+-x64-linux.tar.gz from the NCBI blast website was unzipped to 
# folder /anypath/
# Mus_musculus_.GRCm38.70.cnda.all.fa.gz from the ENSEMBL website was downloaded 
# to /anypath/db/ , extracted, and used to create a blast database.
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
  read.delim("./data/annotation.files/AffyProbe-to-Ensembl-transcript-ID.txt",header=FALSE)
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

# MGI probe sets
CustomPS[["probe"]][["MGI"]] = TransitiveMapping( CustomPS$Ensembl$probe, MsAnn$Ensembl$MGI )
CustomPS[["MGI"]][["probe"]] = inverseList( CustomPS[["probe"]][["MGI"]] )


save(CustomPS,file="./data/annotation.RData/CustomPS.RData")



