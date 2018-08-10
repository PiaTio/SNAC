###Sparse Network And Component analysis
##Integrating two information sources using component and network analysis.  
##
##data: Data file with subjects as rows and variables as columns. Variables should be grouped per data block
##b1: The number of variables in block 1
##b2: The number of variables in block 2
##R: The number of components
##lasso: lambda for graphical lasso. Default is set at 0.002
##NStart: The number of multistarts for the cross-validation algorithm. Default is set to 1. 
##method: "datablock" or "component". These are two options with respect to the grouping of the loadings as used in the Group Lasso penalty. 
##If method="component", the block-grouping of the coefficients is applied per component separately. If method = "datablock", the grouping is applied on the 
##concatenated data block, with loadings of all components together. If method is missing, then the "component" method is used by default.
##nfolds: The number of folds. If missing, then 10 fold cross-validation will be performed.
##MaxIter: Maximum number of iterations for the cross-validation algorithm. Default is set to 400. 


SNAC(data, b1, b2, R, lasso = 0.002, NStart = 1, method, nfolds = 10, MaxIter = 400){
  if(is.numeric(b1)==FALSE){
    stop("b1 is not a numeric entry.")
  }
  
  if(is.numeric(b2)==FALSE){
    stop("b2 is not a numeric entry.")
  }
  
  Jk <- c(b1, b2)
  library(RegularizedSCA)
  
  #Find lambda values for lasso and grouplasso using cross-validation
  maxSeq <- maxLGlasso(DATA = data, Jk = Jk, R = R)
  
  para <- cv_sparseSCA(DATA = data, Jk = Jk, R = R, LassoSequence = seq(0.0001, maxLGasso$Lasso, length.out = 50), 
               GLassoSequence = seq(0.00001, maxLGasso$Glasso, length.out = 50), NRSTARTS = NStart, nfolds = nfolds, MaxIter = MaxIter)
  
  #Run SCA
  results <- sparseSCA(DATA = data, Jk = Jk, R = R, 
            LASSO = para$RecommendedLambda[1], 
            GROUPLASSO = para$RecommendedLambda[2], 
            NRSTART = NStart, method = method, MaxIter = MaxIter)
  
  #Undo shrinkage done on the P and T matrices 
  final <- undoShrinkage(DATA = data, R = R, Phat = results$Pmatrix)
  Pmatrix <- final$Pmatrix
  Tmatrix <- final$Tmatrix
  
  #Find common components
  common = NULL
  for (c in 1:ncol(final$Pmatrix)){ #look per component whether it is a common or distinctive one
    s = 0
    d = NULL
    for (j in 1:length(Jk)){
      d[j] = sum(abs(final$Pmatrix[(s+1):(s+Jk[j]),c]))
      if(d[j]!=0) {d[j] = 1}
      s = sum(Jk[1:j])
    }
    if (sum(d)>= 2) {common = c(common, c)} #common component is defined as having variables of at least two blocks
  }
  TmatCommon = final$Tmatrix[,common]
  PmatCommon = final$Pmatrix[,common]
  
  #Create data set based on common components
  Xcom = TmatCommon%*%t(PmatCommon)
  
  #Graphical lasso
  SigmaINVC_hat = glasso(s = cov(Xcom), rho = lasso, penalize.diagonal = FALSE)$wi
  
  return(SigmaINVC_hat)
}