tauestimate=function(X,Fy,d,aux,inic)
{
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(Fy)[2]
  c1=aux$c1
  c2=aux$c2
  k1=aux$k1
  k2=aux$k2

  #####################################################
  delta0=inic$delta0
  beta0=inic$beta0

  ## Step 1
  aux=X-matrix(c(rep(1,n),Fy),n,(r+1))%*%t(beta0)
  dy=Re(sqrt(diag((aux)%*%solve(delta0)%*%t(aux))))
  s=MscaleNR(dy,c1,k1)$s
  dys=dy/s

  GG=NULL
  G0=phi1(delta0,c2,k2,s,dy,p)  # initial value of the function to be minimized
  GG[1]=G0

  w=weights2(dys,s,p,c1,c2)
  W=mean(w)

  ### Step 2
  W.matrix=matrix(rep(0,n^2),n,n)
  diag(W.matrix)=w
  unos=rep(1,n)
  Xmean=t(X)%*%W.matrix%*%unos/(n*W)
  XX=X-matrix(rep(Xmean,n),n,p,byrow=T)

  Fmean=t(Fy)%*%W.matrix%*%unos/(n*W)
  FF=Fy-matrix(rep(Fmean,n),n,r,byrow=T)

  cov_ff=t(FF)%*%W.matrix%*%FF/(n*W)
  cov_xf=t(XX)%*%W.matrix%*%FF/(n*W)
  pii=Re(cov_xf%*%solve(cov_ff)%*%t(cov_xf))


  ### Step 3
  mues=(t(X)%*%W.matrix%*%unos-beta0[,2:(r+1)]%*%t(Fy)%*%W.matrix%*%unos)/(n*W)

  ## Step 4
  aux2=eigen(sqrtm2(solve(delta0))%*%pii%*%sqrtm2(solve(delta0)))

  gamaes=Re(sqrtm2(delta0)%*%aux2$vectors)[,1:d]

  betaes=Re(solve(t(gamaes)%*%solve(delta0)%*%gamaes)%*%t(gamaes)%*%solve(delta0)%*%t(XX)%*%W.matrix%*%FF%*%solve(t(FF)%*%W.matrix%*%FF))

  ### we improve the weights before estimating delta
  aux=X-matrix(c(rep(1,n),Fy),n,(r+1))%*%t(cbind(mues,gamaes%*%betaes))

  delta.sc=t(aux)%*%W.matrix%*%aux/n

  ## Step 5

  dy=Re(sqrt(diag((aux)%*%solve(delta.sc)%*%t(aux))))

  s=MscaleNR(dy,c1,k1)$s
  dys=dy/s
  cte=mean(rhoc(dys,c2))*s^2/k2
  deltaes=delta.sc*cte

  G1=phi1(deltaes,c2,k2,s,dy,p)  # initial value of the function to be minimized
  GG[2]=G1
  k=2
  while ((abs(log(G0)-log(G1))>1e-5)& k<=50)
  {
    k=k+1
    G0=G1

    ## Step 1
    aux=X-matrix(c(rep(1,n),Fy),n,(r+1))%*%t(cbind(mues,gamaes%*%betaes))
    dy=Re(sqrt(diag((aux)%*%solve(deltaes)%*%t(aux))))
    s=MscaleNR(dy,c1,k1)$s
    dys=dy/s

        w=weights2(dys,s,p,c1,c2)
    W=mean(w)

    ### Step 2
    W.matrix=matrix(rep(0,n^2),n,n)
    diag(W.matrix)=w
    unos=rep(1,n)
    Xmean=t(X)%*%W.matrix%*%unos/(n*W)
    XX=X-matrix(rep(Xmean,n),n,p,byrow=T)

    Fmean=t(Fy)%*%W.matrix%*%unos/(n*W)
    FF=Fy-matrix(rep(Fmean,n),n,r,byrow=T)

    cov_ff=t(FF)%*%W.matrix%*%FF/(n*W)
    cov_xf=t(XX)%*%W.matrix%*%FF/(n*W)
    pii=Re(cov_xf%*%solve(cov_ff)%*%t(cov_xf))

    ### Step 3
    mues=(t(X)%*%W.matrix%*%unos-gamaes%*%betaes%*%t(Fy)%*%W.matrix%*%unos)/(n*W)

    ## Step 4
    aux2=eigen(sqrtm2(solve(deltaes))%*%pii%*%sqrtm2(solve(deltaes)))

    gamaes=Re(sqrtm2(deltaes)%*%aux2$vectors)[,1:d]

    betaes=Re(solve(t(gamaes)%*%solve(deltaes)%*%gamaes)%*%t(gamaes)%*%solve(deltaes)%*%t(XX)%*%W.matrix%*%FF%*%solve(t(FF)%*%W.matrix%*%FF))


    aux=X-matrix(c(rep(1,n),Fy),n,(r+1))%*%t(cbind(mues,gamaes%*%betaes))

    delta.sc=t(aux)%*%W.matrix%*%aux/n

    ## Step 5

    dy=Re(sqrt(diag((aux)%*%solve(delta.sc)%*%t(aux))))

    s=MscaleNR(dy,c1,k1)$s
    dys=dy/s
    cte=mean(rhoc(dys,c2))*s^2/k2
    deltaes=delta.sc*cte
    dy=sqrt(diag((aux)%*%solve(deltaes)%*%t(aux)))
    G1=phi1(deltaes,c2,k2,s,dy,p)  # initial value of the function to be minimized
    GG[k]=G1
      }
  alerta=(k>50)
  if (alerta) message ("maximun iterations are reached")
  #-------
  list(mu=mues,beta=betaes,gamma=gamaes,delta=deltaes)
}



#################################################################
#' computes weights for the algorithm
#'
#' @param dys vector of distances
#' @param s real positive number, the s-cale of distances
#' @param p dimension of ambient space
#' @param c1 constant for rho_1 function
#' @param c2 constant for rho_2 function
#' @NoRd
weights2=function(dys,s,p,c1,c2)
{
  a=2*rhoc(dys,c2)-psic(dys,c2)*dys
  b=psic(dys,c1)*dys
  w1=p*mean(a)/(2*mean(rhoc(dys,c2))*s^2*mean(b))
  w2=p/(2*mean(rhoc(dys,c2))*s^2)
  w=w1*psic(dys,c1)/dys+w2*psic(dys,c2)/dys
  w
  ##----
}


#' computes phi1 function, internally used to compute the weights for the algorithm
#'
#' @param delta positive definite p x p matrix, the estimated covariance of errors
#' @param c2 constant for rho_2 function
#' @param k2 constant for rho_2 function
#' @param ese real positive number, the s-cale of distances
#' @param dis vector of distances
#' @param p dimension of ambient space
#' @NoRd
phi1=function(delta,c2,k2,ese,dis,p)
{
  det(delta)*(mean(rhoc(dis/ese,c2))*ese^2)^p
}


#' computes the principal square root of the matrix A
#'
#' @param A positive (semi) definite p x p matrix
#' @NoRd
sqrtm2=function(A)
{
  # X = sqrtm2(A) is the principal square root of the matrix A, i.e. X\%*\%X = A
  descoA<-eigen(A,symmetric= TRUE)
  descoA$vectors%*%diag(sqrt(descoA$values))%*%t(descoA$vectors)
  #------------
}

