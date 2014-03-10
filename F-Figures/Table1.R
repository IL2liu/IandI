# Read and analyzer the CyberT files

# set up the file names
comparisons = list( c("A2","Sham2"), c("B2","Sham2"), c("AB2","Sham2"),
                    c("A6","Sham6"), c("B6","Sham6"),c("AB6","Sham6"),
                    c("A16","Sham16"), c("B16","Sham16") )
comp.names = sapply(comparisons,function(x) paste(x,collapse="-"))

# get the p-values and fold changes
load("./data/RData/cyberT.results.RData")
pvals = cyberT.results$pval[,comp.names]
fold = cyberT.results$logFC[,comp.names]

# Get probe sets with q<0.01 for TcdA and TcdB
qvals = apply(pvals,2,function(x) p.adjust(x,"BH"))
a.probesets = names( which(qvals[,c("A2-Sham2")]<0.01) )
b.probesets = names( which(qvals[,c("B2-Sham2")]<0.01) )

# Convert Affy ID to gene symbol
load("./data/annotation.RData/MsAnn.RData")
a.genes = sapply( MsAnn$Affy$Symbol[a.probesets], function(x) x[1] )
b.genes = sapply( MsAnn$Affy$Symbol[b.probesets], function(x) x[1] )

# Some probe sets don't map to any genes. Their symbol is NULL
# and they are removed by the unlist function.
a.genes = unlist(a.genes)
b.genes = unlist(b.genes)

# Now save the probe sets which mapped to a gene
a.probesets = names(a.genes)
b.probesets = names(b.genes)

# and make tables with these probe sets
a.fold = fold[a.probesets,]
rownames(a.fold) = a.genes
b.fold = fold[b.probesets,]
rownames(b.fold) = b.genes

# sort by gene expression at 2 hours
a.fold = a.fold[order(-a.fold[,1]),]
b.fold = b.fold[order(-b.fold[,2]),]

# combine the two tables
tab = rbind(a.fold,b.fold)

# "unlog" the fold change. This is Table 1
table1 = round( sign(tab)*2^(abs(tab)), 1)
table1


