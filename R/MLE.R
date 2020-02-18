MLE=function(X,Fy,d)

{
  if (is.vector(X)) n=length(X) else n=dim(X)[1]
  if (is.vector(X)) p=1 else p=dim(X)[2]
  Xmean=apply(X,2,mean)
  XX=t(t(X)-rep(Xmean,n))
  Fy=as.matrix(Fy)


  Fmean=apply(Fy,2,mean)

  FF=t(t(Fy)-rep(Fmean,n))

  Sigmahat=t(XX)%*%XX/n
  Sigmafit=t(XX)%*%FF%*%solve(t(FF)%*%FF)%*%t(FF)%*%XX/n
  Sigmares=Sigmahat-Sigmafit

  A=solve(Sigmares)
  V=eigen(sqrtm2(A)%*%Sigmafit%*%sqrtm2(A))$vectors
  autoval=eigen(sqrtm2(A)%*%Sigmafit%*%sqrtm2(A))$values
  autoval[1:d]=0
  K=diag(autoval)

  Delta=Sigmares+sqrtm2(Sigmares)%*%V%*%K%*%t(V)%*%sqrtm2(Sigmares)
  B=solve(Delta)
  #B=Delta
  aux1=eigen(sqrtm2(B)%*%Sigmafit%*%sqrtm2(B))$vectors
  Gama=sqrtm2(Delta)%*%aux1[,1:d]

  beta=solve(t(Gama)%*%B%*%Gama)%*%t(Gama)%*%B%*%t(XX)%*%FF%*%solve(t(FF)%*%FF)
  mu=Xmean-Gama%*%beta%*%Fmean

  list(mu=mu,beta=beta,gamma=Gama,delta=Delta)
  #----
}

