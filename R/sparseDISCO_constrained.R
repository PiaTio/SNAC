sparseDISCO_constrained<-function(Xc,Pinit)
{
#sparseDISCO results in sparse common and distinctive simultaneous components
#INPUT:
#Xc:			concatenated data
#JkSize:		number of variables in each data block
#Rcommm_dist:	number of components of each kind (common, distinctive for each block; in this order!)
#LL:			value of lasso tuning parameter
#LG:			value of group lasso tuning parameter (LARGE?)
#LE:			value of elitist lasso tuning parameter

#initialization
d<-dim(Pinit)
J<-d[1]
I<-nrow(Xc)
R<-d[2]
P<-matrix(rnorm(J*R),nrow=J,ncol=R,byrow=TRUE)
P[Pinit==0]<-0
EPS<-1e-16#rounding off bound

#alternating routine
conv<-0
iter<-0
#P<-Pinit
while (conv==0)
	{
	#conditional estimation of T given P: closed form
      A<-t(Xc%*%P)
      svdA<-svd(A)
      T<-svdA$v%*%t(svdA$u)
	
	#check non-increasing loss
	DEV<-Xc-(T%*%t(P))
	DEVsq<-DEV^2
	fit<-sum(apply(DEVsq,2,sum))
	Loss<-fit

cat("Update T Loss: ",Loss,sep="\n")

	#conditional estimation of P given T with some P fixed to zero
	for (r in 1:R)
	{
		CPr=t(Xc)%*%T[,r]
		P[,r]=CPr
	}               
	P[Pinit==0]<-0
             	
	#check non-increasing loss
	DEV<-Xc-(T%*%t(P))
	DEVsq<-DEV^2
	fit<-sum(apply(DEVsq,2,sum))
	Loss<-fit

cat("Update P Loss: ",Loss, sep="\n")

iter=iter+1
if(MAXITER==iter) conv<-1
}
Loss<-Loss/sum(rowSums(Xc^2))
  result <- list(T=T,P=P,L=Loss)
  return(result)
}#end of function