angle_AB=function(A,B)
{ #compute angles between subspaces where the columns of A and B are the respective bases
 
  A=as.matrix(A)
  B=as.matrix(B)
  d=min(ncol(A),ncol(B))
  P.A=(A)%*%solve(t(A)%*%A)%*%t(A)
  P.B=(B)%*%solve(t(B)%*%B)%*%t(B)
  aa<-(eigen(t(P.A)%*%P.B%*%P.A,symmetric= TRUE)$values)[1:d]
  aa[(aa-1)>10^(-18)]<-rep(1,sum(((aa-1)>10^(-18))))
  acos(sqrt(aa))*180/pi
  ##------
}