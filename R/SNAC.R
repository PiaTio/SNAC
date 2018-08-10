###Sparse Network And Component analysis
##Integrating multi-source informatie using component and network analysis 


SNAC(data, b1, b2, R, stand = TRUE, lasso = 0.002, NStart, method){
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
               GLassoSequence = seq(0.00001, maxLGasso$Glasso, length.out = 50), NRSTARTS = NStart)
  
  #Run SCA
  results <- sparseSCA(DATA = data, Jk = Jk, R = R, 
            LASSO = para$RecommendedLambda[1], 
            GROUPLASSO = para$RecommendedLambda[2], 
            NRSTART = NStart, method = method)
  
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