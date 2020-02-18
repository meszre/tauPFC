tauRRR=function(yy,XX,d,efficiency=0.90)
{
  p=dim(yy)[2]
  r=dim(XX)[2]

  alerta=(d>min(p,r))
  if (alerta) stop ("rank should be smaller or equal to min(dim(yy)[2],dim(XX)[2])")

  aux1<-kappa_and_c(p,efficiency = 0.85)
  # choose 0.85 efficiency to control bias of the initial estimator
  # keeping efficiency not to low
  primero<-initial(X = yy, Fy = XX, aux = aux1, efficiency = 0.85)
  aux2<-kappa_and_c(p,efficiency)
  segundo<-tauestimate(X = yy, Fy = XX, d, aux = aux2, inic = primero)
  # in the notation of Izenman (2008), under the RRR model,
  # the robust estimate of the vectorial intercept is
  mu=segundo$mu

  # in the notation of Izenman (2008), a possible A matrix is
  BB=segundo$beta
  # in the notation of Izenman (2008), a possible B matrix is
  AA=segundo$gamma


  # then the reduced-rank regression coefficient matrix
  # of dimension p times dim(XX)[2] with rank = d is
  RRRcoef=AA%*%BB

  # the estimate of the covariance matrix of the errors of the
  # MLM is

  cov.error=segundo$delta

  list(mu=mu,RRRcoef=RRRcoef,AA=AA,BB=BB,cov.error=cov.error)
}
