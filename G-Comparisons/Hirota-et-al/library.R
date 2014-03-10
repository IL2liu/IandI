
################################################################################
#                 Functions for parsing data                                   #
################################################################################

getSerumData = function() {
  raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
  raw = raw[ !is.na(raw$Cxcl1.Serum) &
               !is.na(raw$Cxcl2.Serum) , 
             c("ID","Time","Toxin","Antibody","Death",
               colnames(raw)[grep("Serum",colnames(raw))])]
  raw$Death = factor(as.numeric(raw$Death!=""))
  
  id.vars = c("ID","Toxin","Time","Antibody","Death")
  prot.vars = grep("[.]Serum",colnames(raw),value=TRUE)
  proteins = gsub("[.]Serum","",prot.vars)
  raw = na.omit(raw[,c(id.vars,prot.vars)])
  colnames(raw) = c(id.vars, proteins )
  mraw = melt(raw,id.vars=id.vars)
  avs = cast( mraw, Toxin+Antibody ~ variable, mean )
  lavs = log2(as.matrix.data.frame(avs[,3:12]))
  dlavs = rbind( lavs[1,] - lavs[3,], lavs[2,] - lavs[4,] )
  lfc = cbind( avs[1:2,1:2], dlavs )
  
  lfc2 = lfc[,3:12]
  rownames(lfc2) = c("A.Anticxcl","A.Isotype")
  avs2 = avs[,3:12]
  rownames(avs2) = c("A.Anticxcl","A.Isotype","Sham.Anticxcl","Sham.Isotype")
  
  return(list(raw=avs2,lfc=lfc2))
}

getHirotaData = function() {
  aa = read.csv("./G-Comparisons/Hirota-et-al/Cytokines-colonicTissueLysates.csv")
  bb = aa[,c(2,3,5)]
  rownames(bb) = bb[,1]
  bb[ bb == "BDL" ] = NA
  cc = bb[,-1]
  cc[,1] = as.numeric(as.character(cc[,1]))
  cc[,2] = as.numeric(as.character(cc[,2]))
  dd = t(cc)
  lfc = log2(dd)[2,] - log2(dd)[1,]
  lfc2 = matrix(lfc,nrow=1)
  colnames(lfc2) = names(lfc)
  rownames(lfc2) = "AB"
  
  return(list(raw=dd,lfc=lfc2))
}

getCellData = function() {
  
  # Now compare this to our protein data from cells
  raw = read.csv("./PhysiologyMeasurements.csv",sep="\t")
  id.vars = c("ID","Toxin","Time","Antibody","Death")
  prot.vars = grep("[.]cell",colnames(raw),value=TRUE)
  proteins = gsub("[.]cell","",prot.vars)
  raw = na.omit(raw[,c(id.vars,prot.vars)])
  raw$Death = factor(as.numeric(raw$Death!=""))
  colnames(raw) = c(id.vars, proteins )
  mraw = melt(raw,id.vars=id.vars)
  
  # average replicates and get fold change
  avs = cast( mraw, Toxin ~ variable, mean )
  rownames(avs) = avs[,1]
  log.avs = log2(as.matrix.data.frame(avs[,-1]))
  prot.lfc = sweep(log.avs,2,log.avs[3,])[1:2,]
  
  avs2 = avs[,-1]
  rownames(avs2) = avs[,1]
  
  return(list(raw=avs2,lfc=prot.lfc,data=raw))
}


getGeneData = function(symbols) {
  load("./DiffProbeSets/data/RData/cyberT.results.RData")
  load("./data/annotation.RData/MsAnn.RData")
  load("./data/annotation.RData/OrthoAnn.RData")
  source("./C-Annotations/TransitiveMapping.R")
  map = TransitiveMapping( MsAnn$MGI$Symbol, MsAnn$MGI$Ensembl )
  ensgs = unlist(map[symbols])
  lfc = cyberT.results$logFC[ensgs,]
  rownames(lfc) = symbols
  return(lfc)
}

corplot = function(data) {
  heatmap.2(cor(data),cellnote=round(cor(data),3),
            col=bluered,symbreaks=TRUE,notecol="black",trace="none")
}
