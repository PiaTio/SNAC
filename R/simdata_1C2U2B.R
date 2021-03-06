###Sparse Network And Component analysis
##Generating 2-block gaussian distributed data with underlying 3 component structure (1 common, 2 unique)
##
##m: The desired number of participants.
##b1: The number of variables in block 1.
##b2: The number of variables in block 2.
##sp1: Desired sparseness in component structure for first data block. Default set to 0.5.
##sp2: Direced sparseness in component structure for second data block. Default set to 0.8.
##S: Indicates whether the distinctives components are dominanting the overall component structure (2) or not (1). Default is set to not (1). 

simdata_1C2U2B <- function(m, b1, b2, sp1 = 0.5, sp2 = 0.8, S = 1){
  
  #Number of participants
  m <- m
  #Block sizes
  JkSize<-c(b1, b2)
  #Number of variables
  n <- sum(JkSize)
  #Number of common and unique for each block components
  Rcomm_dist<-c(1,2)
  #Number of componenten
  R <- sum(Rcomm_dist)
  
  #### Singular value decomposition ####
  #Basic data input
  X <- matrix(rnorm(m*n),nrow=m)
  #Standardise
  Z <- scale(X)
  #SVD
  svdX <- svd(X)
  #Common subjects --> same U
  U <- svdX$u[,1:R]
  #Common subjects --> different V
  Vcomm <- svdX$v[,1:Rcomm_dist[1]]
  #How sparse is common component
  sp <- 0.5 
  Zc <- sample(c(1:Rcomm_dist[1]*n),sp*n)
  Vcomm[Zc] <- 0
  #How sparse are the unique components
  Vdist <- svdX$v[,(Rcomm_dist[1]+1):R]
  #Create first block
  Vdist[1:JkSize[1],1] <- 0
  #How spare is unique component 1
  sp1 <- sp1
  Zu1 <- sample(c((1+JkSize[1]):(JkSize[1] +JkSize[2])),sp1*JkSize[1])
  Vdist[Zu1, 1] <- 0
  #Create second block
  Vdist[(1+JkSize[1]):(JkSize[1] +JkSize[2]),2] <- 0
  #How sparse is unqie component 2
  sp2 <- sp2
  Zu2 <- sample(1:JkSize[1], sp2*JkSize[2])
  Vdist[Zu2, 2] <- 0
  #Create V with 1 common component (mixed with 0), one unique component for block 2 and a
  #one unique component for block 1
  Vstar <- cbind(Vcomm,Vdist)
  
  #create S (singular values)
  s1 <- 1
  s2 <- 1
  s3 <- 1
  Sstar <- matrix(data = 0, nrow = R, ncol = R)
  diag(Sstar) <- c(s1, s2, s3)
  #Adjust V based on S
  a <- matrix(nrow = n, ncol = R)
  for (i in 1:R){
    a[,i] <- (Vstar[,i]^2)*Sstar[i,i]^2
  }
  A <- rowSums(a)
  
  #Adjust V to make sure that correlation Ris between -1 and 1, and imposing correlation structure
  # Here block 1 r = 0.8 and block 2 r = 0.3
  Vdstar <- matrix(0, sum(JkSize), R)
  Vdstar[1:100,] <- sqrt(0.8/A[1:100])*Vstar[1:100,] 
  Vdstar[101:200,] <- sqrt(0.3/A[101:200])*Vstar[101:200,]  
  Vdstar[Vdstar=="NaN"]<-0 
  #True covariance structure
  cov <- Vdstar%*%(Sstar^2)%*%t(Vdstar) 
  diag(cov)<-1
  library(MASS)
  inv <- solve(cov)
  ginv <- ginv(cov)
  #Calculate population covariance matrix for common component
  cov23 <- Vdstar[,2:3]%*%(Sstar[2:3,2:3]^2)%*%t(Vdstar[,2:3])
  cov1 <- cov-cov23
  inv1 <- solve(cov1)
  ginv1 <- ginv(cov)
  
  #Data
  X <- mvrnorm2(m, rep(0, n), cov, empirical = TRUE)
  
  return(list(X, cov, inv, ginv, cov1, inv1, ginv1, Vdstar, Sstar, Vstar))
  
}


