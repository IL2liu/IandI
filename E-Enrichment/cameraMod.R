library(limma)

################################################################################
#                            modified CAMERA                                   #
################################################################################

# CAMERA modified to use a different bayesian test statistic. 
# The prior degrees of freedom must be given.
# see ?camera for more information (limma package)
# tStat is a vector of the test statistics
# d0 is the prior degrees of freedom
cameraMod = function( indices, y, design, contrast=ncol(design), 
                        tStat, d0 ) {
  
  ## Getting the necessary variables
  G = length(tStat)
  n = nrow(design)
  p = ncol(design)
  df.residual = n-p
  df.camera = min(df.residual, G-2)
  contrast = ncol(design)
  j = c((1:p)[-contrast], contrast)
  
  # getting the log fold change
  QR <- qr(design)
  effects <- qr.qty(QR, t(y))
  unscaledt <- effects[p, ]
  
  # Getting the variance of each gene
  U = effects[-(1:p), , drop = FALSE]
  sigma2 = colMeans(U^2)
  U = t(U)/sqrt(sigma2)
  A = NULL
  
  # Total degrees of freedom for the t-distribution of the t-statistic
  # used to make scores from the modified t-statistics
  df.total = min(df.residual + d0, G * df.residual)
  Stat = zscoreT(tStat, df = df.total)
  meanStat <- mean(Stat)
  varStat <- var(Stat)
  
  # prepare results table
  nsets <- length(indices)
  tab <- matrix(0, nsets, 5)
  rownames(tab) <- names(indices)
  colnames(tab) <- c("NGenes", "Correlation", "Down", "Up", 
                     "TwoSided")
  
  for (i in 1:nsets) {
    
    # prepare variables for set level test
    index <- indices[[i]]
    StatInSet <- Stat[index]
    m <- length(StatInSet)
    m2 <- G - m
    
    # calculate the variance inflation factor
    Uset <- U[index, , drop = FALSE]
    vif <- m * mean(colMeans(Uset)^2)
    correlation <- (vif - 1)/(m - 1)
    
    # the variance inflated t-test for enrichment
    meanStatInSet <- mean(StatInSet)
    delta <- G/m2 * (meanStatInSet - meanStat)
    varStatPooled <- ((G - 1) * varStat - delta^2 * m * 
                        m2/G)/(G - 2)
    two.sample.t <- delta/sqrt(varStatPooled * (vif/m + 
                                                  1/m2))
    
    # fill in the results table
    tab[i, 1] <- m
    tab[i, 2] <- correlation
    tab[i, 3] <- pt(two.sample.t, df = df.camera)
    tab[i, 4] <- pt(two.sample.t, df = df.camera, lower.tail = FALSE)
    
    
  }
  tab[, 5] <- 2 * pmin(tab[, 3], tab[, 4])
  return(tab)
}








