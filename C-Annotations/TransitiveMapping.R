library(lattice)
library(Matrix)

list2matrix = function( list1 ) {
  out = as.matrix( cbind( rep(names(list1),times=sapply(list1,length)),
                   unlist(list1)) )
  dimnames(out)=NULL
  return(out)
}

inverseList = function( inlist ) {
  outlist = split( rep(names(inlist),times=sapply(inlist,length)),  unlist(inlist) )
  return(outlist)
}


# This function connects two separate maps of identifiers that have a common
# identifier between them. If, for instance, one matrix maps identifiers x to a 
# and another maps identifers x to b, then this function maps identifers a
# to b. This mapping is performed efficiently by using matrix multiplication.
TransitiveMapping = function(IDx.to.IDa,IDx.to.IDb) {
  
  if( typeof(IDx.to.IDa)=="list" && typeof(IDx.to.IDb)=="list"){
    IDx.to.IDa = list2matrix(IDx.to.IDa)
    IDx.to.IDb = list2matrix(IDx.to.IDb)
    listflag = "yes"
  } else {
    listflag = "no"
  }
  
  # map x to a
  IDx1 = IDx.to.IDa[,1]
  IDa = IDx.to.IDa[,2]
  
  # map x to b
  IDx2 = IDx.to.IDb[,1]
  IDb = IDx.to.IDb[,2]
  
  # get the union of all the x IDs
  IDx = unique(c(IDx1,IDx2))
  
  # give each ID in x a number
  IDcx = 1:length(IDx)
  names(IDcx) = as.character(IDx)
  
  # give each ID in "a" a number
  IDca = 1:length(unique(IDa))
  names(IDca) = as.character(unique(IDa))
  
  # give each ID in "b" a number
  IDcb = 1:length(unique(IDb))
  names(IDcb) = as.character(unique(IDb))
  
  # make a matrix representing map x to a
  top1 = as.numeric(IDcx[as.character(IDx1)])
  left1 = IDca[as.character(IDa)]
  x2a = sparseMatrix(i=left1,j=top1,dims=c(length(IDa),length(IDcx)) )
  
  # make a matrix representing map x to b
  top2 = as.numeric(IDcx[as.character(IDx2)])
  left2 = IDcb[as.character(IDb)]
  x2b = sparseMatrix(i=left2,j=top2,dims=c(length(IDb),length(IDcx)) )
  
  # matrix multiplication
  final = x2a %*% t(x2b)
  
  # Go backwards from numbers just created to IDs
  left = summary(final)$i
  right = summary(final)$j
  
  nleft = names(IDca)[left]
  nright = names(IDcb)[right]
  aa = cbind(nleft,nright)
  
  if( listflag=="yes" ) {
    return(split(aa[,2],aa[,1]))
  } else {
    return(aa)
  }
  
}


