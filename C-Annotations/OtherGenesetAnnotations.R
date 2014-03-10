# Parse the file to link Human ENTREZ IDs to REACTOME
# gene assocations

# The container
GSAnn = list()

# Parsing the file
aa = as.matrix(read.delim(
  "./data/annotation.files/REACTOME-to-Human-Entrez.gmt",sep="\n",header=F))
foo = function(x) strsplit(x,"\t")
bb = sapply(aa,foo)
names(bb) = NULL
set.names = lapply(bb,function(x) x[1])
sets = lapply(bb,function(x) x[3:length(x)])
names(sets) = set.names

# Saving to the "containers"
GSAnn[["REACTOME"]][["HumanENTREZ"]] = sets
entrez.reactome = split( rep(names(sets),times=sapply(sets,length)),
                         unlist(sets) )
GSAnn[["HumanENTREZ"]][["REACTOME"]] = entrez.reactome

# file
save(GSAnn, file="./data/annotation.RData/GSAnn.RData")



