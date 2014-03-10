# A script that runs cyberT for all the different two class comparisons
load("./DiffProbeSets/data/RData/esets.RData")
source("./D-DiffExpression/cybertWrapper.R")

# Arrays were preprocessed in groups designated by time point (e.g. 2h, 6h or 16h)
# Here combine them only for the sake of simpler code
eset = new("ExpressionSet")
e1 = esets[["2hr"]]
e2 = esets[["6hr"]]
e3 = esets[["16hr"]]
exprs(eset) = cbind(exprs(e1),exprs(e2),exprs(e3))
pData(eset) = rbind(pData(e1),pData(e2),pData(e3))

comparisons = list( c("A2","Sham2"), c("B2","Sham2"), c("AB2","Sham2"),
                    c("A6","Sham6"), c("B6","Sham6"), c("AB6","Sham6"),
                    c("A16","Sham16"), c("B16","Sham16") )

# Run cyberT and save the results
cyberT.results = cybertWrapper( eset, comparisons )
save(cyberT.results,file="./DiffProbeSets/data/RData/cyberT.results.RData")
