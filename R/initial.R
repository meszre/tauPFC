initial=function(X,Fy,aux,efficiency=0.85)
{
  p=dim(X)[2]
  r=dim(Fy)[2]
  n=dim(X)[1]
  c1=aux$c1
  c2=aux$c2
  k1=aux$k1
  k2=aux$k2
  options(warn = -1)
  #library(robust)
  #library(rrcov)
# if (!requireNamespace("robust", quietly = TRUE)) {
#  stop("Package \"robust\" needed for this function to work. Please install it.",
#         call. = FALSE)
#  }

  BETA=matrix(nrow=p,ncol=r+1)
  for (kkk in 1:p){
    BETA[kkk,]=robust::lmRob(X[,kkk]~Fy,control = robust::lmRob.control(efficiency = efficiency))$coefficients  ## library(robust)
  }
  res=X-cbind(rep(1,n),Fy)%*%t(BETA)
  A=rrcov::CovSest(res,method="bisquare")  #library(rrcov)
  delta0.aux=rrcov::getCov(A)

  aux2=X-matrix(c(rep(1,n),Fy),n,(r+1))%*%t(BETA)
  dy=sqrt(diag((aux2)%*%solve(delta0.aux)%*%t(aux2)))
  s= MscaleNR(dy,c1,k1)$s

  dys=dy/s
  cte=mean(rhoc(dys,c2))*s^2/k2
  delta0=delta0.aux*cte
  list(beta0=BETA, delta0=delta0)

}


#' computes the difference between the mean of rhoc values and 0.5
#'
#' @param c1 constant for rho_1 function
#' @param u vector of distances (positive numbers)
#' @param s real positive number, the s-cale of distances
#' @NoRd
h=function(s,u,c1){mean(rhoc(u/s,c1))-.5 }

#' computes Mscale using Newton Raphson
#'
#' @param c1 constant for rho_1 function
#' @param k1 constant for rho_1 function
#' @param u vector of distances (positive numbers)
#' @NoRd
MscaleNR<-function(u,c1,k1){
  sn=stats::median(abs(u))/.6745
  while (h(sn,abs(u),c1)>0){sn=1.5*sn}
  i<-0
  err <- 1
  while  (( i < 1000 ) & (err > 1e-5)) {
    var1=  u/sn;
    A=mean(rhoc(var1,c1));
    B=mean(psic(var1,c1)*var1);
    factorAB=(A -B -0.5)/(2*A-B-1);
    snmas1 = sn*factorAB;
    err <- abs(snmas1/sn - 1)
    sn=abs(snmas1)
    i  <- i+1
  }
  return(list(s=snmas1,numeroDeIteraciones=i))
}
