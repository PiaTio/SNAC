sparseDISCO_knownnrzeros<-function(Xc,JkSize,Rcomm_dist,LLd,LLc,MAXITER,INIT,Nr0)
{
  #sparseDISCO results in sparse common and distinctive simultaneous components
  #INPUT:
  #Xc:			concatenated data
  #JkSize:		number of variables in each data block
  #Rcommm_dist:	number of components of each kind (common, distinctive for each block; in this order!)
  #LLd:			value of lasso tuning parameter for distinctive components
  #LLc:			value of lasso tuning parameter for common components
  #INIT:		to be used in case of known zero coefficients: matrix of size nr.of vars x nr. of components; zero values remain zero
  #Nr0:			desired number of zero coefficients
  
  #initialization
  J<-sum(JkSize)
  I<-nrow(Xc)
  K<-length(JkSize)
  R<-sum(Rcomm_dist)
  svdX<-svd(Xc,nv=R)
  P<-svdX$v
  
  #P<-matrix(rnorm(J*R),nrow=J,ncol=R,byrow=TRUE)
  if(prod(dim(INIT))!=0)
  {
    P[INIT==0]<-0
  }
  #weight penalty for each coefficient: randomized lasso
  #B<-matrix(runif(J*R),nrow=J,ncol=R,byrow=TRUE)
  #B[B<0.2]<-0.2
  B<-matrix(c(rep(1,J*R)),nrow=J,ncol=R,byrow=TRUE)
  EPS<-1e-12#rounding off bound
  
  penL<-0
  
  #COMMON part
  L<-0
  if(Rcomm_dist[1]>0){
    rcd<-Rcomm_dist[1]
    for (k in 1:K)
    {
      Jk<-JkSize[k]
      subP<-matrix(P[(L+1):(L+Jk),1:rcd],ncol=rcd)
      subB<-matrix(B[(L+1):(L+Jk),1:rcd],ncol=rcd)
      penL<-penL+sum(colSums(abs(subB*subP)))*LLc
      L<-L+Jk
    }
  }
  
  #DISTINCTIVE part
  PGmat<-c()
  L<-0
  if(Rcomm_dist[2]>0){
    rcd<-c(seq(Rcomm_dist[1]+1,Rcomm_dist[2]+Rcomm_dist[1]))
    for (k in 1:K)
    {
      Jk<-JkSize[k]
      subP<-matrix(P[(L+1):(L+Jk),rcd],ncol=Rcomm_dist[2])
      subB<-matrix(B[(L+1):(L+Jk),rcd],ncol=Rcomm_dist[2])
      penL<-penL+sum(colSums(abs(subB*subP)))*LLd
      L<-L+Jk
    }
  }
  
  #alternating routine
  conv<-0
  iter<-0
  while (conv==0)
  {
    #conditional estimation of T given P: closed form
    A<-(Xc%*%P)
    svdA<-svd(A)
    T<-svdA$u%*%t(svdA$v)
    
    #check non-increasing loss
    DEV<-Xc-(T%*%t(P))
    DEVsq<-DEV^2
    fit<-sum(apply(DEVsq,2,sum))
    Loss<-fit+penL
    
    cat("Update T Loss: ",Loss,sep="\n")
    
    #conditional estimation of P given T: majorization and thresholding (lasso)
    Pold<-P
    for (r in 1:R)
    {
      CPr=t(Xc)%*%T[,r]
      ifelse(r>Rcomm_dist[1],LL<-LLd,LL<-LLc)
      pols=sign(CPr)*(abs(CPr)-LL*B[,r]/2)
      p=(pols);
      p[abs(CPr)<LL*B[,r]/2]=0
      P[,r]=p
    }               
    P[abs(Pold)<EPS]<-0
    
    penL<-0
    
    
    #COMMON part
    L<-0
    PEmat<-c()
    if(Rcomm_dist[1]>0){
      rcd<-Rcomm_dist[1]
      for (k in 1:K)
      {
        Jk<-JkSize[k]
        subP<-matrix(P[(L+1):(L+Jk),1:rcd],ncol=rcd)
        subB<-matrix(B[(L+1):(L+Jk),1:rcd],ncol=rcd)
        penL<-penL+sum(colSums(abs(subB*subP)))*LLc
        L<-L+Jk
      }
    }
    
    #DISTINCTIVE part
    PGmat<-c()
    L<-0
    if(Rcomm_dist[2]>0){
      rcd<-c(seq(Rcomm_dist[1]+1,Rcomm_dist[2]+Rcomm_dist[1]))
      for (k in 1:K)
      {
        Jk<-JkSize[k]
        subP<-matrix(P[(L+1):(L+Jk),rcd],ncol=Rcomm_dist[2])
        subB<-matrix(B[(L+1):(L+Jk),rcd],ncol=Rcomm_dist[2])
        penL<-penL+sum(colSums(abs(subB*subP)))*LLd
        L<-L+Jk
      }
    }
    
    #check non-increasing loss
    DEV<-Xc-(T%*%t(P))
    DEVsq<-DEV^2
    fit<-sum(apply(DEVsq,2,sum))
    Loss<-fit+penL
    
    cat("Update P Loss: ",Loss, sep="\n")
    
    iter=iter+1
    #stop criteria: maximum nr of iterations 
    if(MAXITER==iter) conv<-1
    #stop criteria: desired number of zeros reached
    if(sum(P==0)>=Nr0) conv<-1
  }
  result <- list(T=T,P=P,L=Loss)
  return(result)
}#end of function