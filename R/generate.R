generate=function(p,n,mutrue,gammatrue,betatrue,sigmatrue){
  ##----
  ## from gamma_0 and beta_0 generate samples according Section 7.1-case d=1
  ##require("MASS")
  y=stats::runif(n,0,4)
  r=2
  fy=matrix(c(y,y^2),n,r)
  X=t(mutrue+gammatrue%*%betatrue%*%t(fy))+sigmatrue*MASS::mvrnorm(n,rep(0,p),diag(p))
  list(Fy=fy,X=X)
  ##----
}
