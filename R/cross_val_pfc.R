cross_val_pfc=function(X,fy,dmax,aux,grafico)
{
# five fold cross validation for pfc model, to estimate d,
# the dimension of the reduction subspace, using both maximum likelihood
# and robust tau estimators to select d. Computes sd of both objective
# functions and gives back the smaller d that satisfies the
# one standard deviation rule from the minimum (see Hastie & Tibshirani)

n=dim(X)[1]
p=dim(X)[2]
r=dim(fy)[2]

c1=aux$c1
c2=aux$c2
k1=aux$k1
k2=aux$k2

# X and fy has the observations by row, dim(X)[1] should be equal to dim(fy)[1]

if ( dim(X)[1] != dim(fy)[1]) stop ("dim(X)[1] should be equal to dim(fy)[1], each row should be an observation")

n_parte=floor(n/5)

if (dmax>p) stop ("dmax should be smaller than dim(X)[2]")


totald=length(0:dmax)

# to save the objective function of ML for the folds
objetivos.mv<-array(0,dim=c(5,totald))


# to save the value of the objective function of ML for the sample
obj.mv<-numeric(totald)
obj.mv.sd<-numeric(totald)

# to save the objective function of tau-estimator for the folds
objetivos.rob<-array(0,dim=c(5,totald))


# to save the value of the objective function for tau for the sample
funcion.FiB<-numeric(totald)
funcion.FiB.sd<-numeric(totald)


particion<-sample(1:n, size=n, replace = FALSE, prob = NULL)


for(j in 1:5)
{
  # (j) is the fold (five folds)
  # we select 4/5 parts of the sample to estimate the parameters
  # this is the training data
  training<-particion[-((n_parte*(j-1)+1):(n_parte*(j-1)+n_parte))]

  # the indices of 1/5 of the sample that is used to validate is kept in "validating"
  # we compute the mahalanobis distance for each one of this 1/5 observations
  # and we compute both objective functiones (max lik and robust) in the validation sample

  validating<-particion[(n_parte*(j-1)+1):(n_parte*(j-1)+n_parte)]  #are the indices of the validation part of the sample

  if (j==5){
    # if n is multiple of 5 this is the same, if not, we redefine training and validating sample to use all the sample
    training<-particion[-((n_parte*(j-1)+1):n)]
    validating<-particion[(n_parte*(j-1)+1):n]  #are the indices of the validation part of the sample
  }
  ##############################
  # first consider d = 0
  ##############################

  d=0

  # we now estimate the parameters
  # maximum likelihood
  mucero<-apply(X[training,],2,mean)
  delta_d_0<-stats::var(X[training,])*(length(training)-1)/length(training)

  # robust
  inic<-initial_d0(X[training,],aux,efficiency=0.85)
  aa<-tauestimate_d0(X[training,],aux,inic)
  mucerorob<-aa$mu
  delta_d_0_rob<-aa$delta

  #now we use them to validate
  #maximum likelihood method
  resiMV<-t(X[validating,])-matrix(rep(mucero,length(validating)),nrow=p, ncol=length(validating))
  dismvfold<-sqrt(diag(t(resiMV)%*%solve(delta_d_0)%*%resiMV))

  objetivos.mv[j,(dmax+1)]=mean(dismvfold^2)+log(det(delta_d_0))

  #robust tau method
  resirob<-t(X[validating,])-matrix(rep(mucerorob,length(validating)),nrow=p, ncol=length(validating))
  distaufold<-sqrt(diag(t(resirob)%*%solve(delta_d_0_rob)%*%resirob))

  sfold=MscaleNR(distaufold,c1,k1)$s

  objetivos.rob[j,(dmax+1)]<-log(det(delta_d_0_rob)) + p*log(mean(rhoc(distaufold/sfold,c2)))+ p*log(sfold^2)


  #############################################
  # second, consider d > 0 (for the same fold)
  #############################################

  # we compute the initial robust estimate

  inic=initial(X[training,],fy[training,],aux,efficiency=0.85)
  for (d in 1:dmax){

      # we now estimate the parameters
      # maximum likelihood
      lme=MLE(X[training,],fy[training,],d)

      # robust tau
      tau_aux=tauestimate(X[training,],fy[training,],d,aux,inic)


      # the indices of 1/5 of the sample that is used to validate is kept in "validating"
      # we compute the mahalanobis distance for each one of this 1/5 observations
      # and we compute both objective functiones (max lik and robust) in the validation sample

      #maximum likelihood method
      resiMV<-t(X[validating,])-matrix(rep(lme$mu,length(validating)),nrow=p, ncol=length(validating))-lme$gamma%*%lme$beta%*%t(fy[validating,])
      dismvfold<-sqrt(diag(t(resiMV)%*%solve(lme$delta)%*%resiMV))
      objetivos.mv[j,d]=mean(dismvfold^2)+log(det(lme$delta))



      #robust tau method
      resirob<-t(X[validating,])-matrix(rep(tau_aux$mu,length(validating)),nrow=p, ncol=length(validating))-tau_aux$gamma%*%tau_aux$beta%*%t(fy[validating,])
      distaufold<-sqrt(diag(t(resirob)%*%solve(tau_aux$delta)%*%resirob))

      sfold=MscaleNR(distaufold,c1,k1)$s
      objetivos.rob[j,d]=log(det(tau_aux$delta)) + p*log(mean(rhoc(distaufold/sfold,c2)))+ p*log(sfold^2)

      # end of the fold over the d's
    }

}


# with the objective functions  computed for all d and for each fold,
# we compute the value of both objective functions for each d


for (d in 1:(dmax+1)){

  #maximum likelihood
  obj.mv[d]<- mean(objetivos.mv[1:5,d])


  #robust
  funcion.FiB[d] = mean(objetivos.rob[1:5,d])

  # sd computations, by tibshirani
  obj.mv.sd[d]<-sd(objetivos.mv[1:5,d])/sqrt(5)
  funcion.FiB.sd[d]<-sd(objetivos.rob[1:5,d])/sqrt(5)

}


# finding minimum and applying one-standard-deviation rule, to select d

  indice.min.mv<-which.min(obj.mv)  # is a number between 1 and 6

# then min(obj.mv)==obj.mv[indice.min.mv]
# in dselec.mv we include all dd that satisfy one-standard-deviation
# rule, then we select the minimum (i.e. the most parsimonious model),
# and store it in dselec.mv


dselec.mv<-indice.min.mv
# the idea is to make dselec.mv a vector of possible values under indice.min.mv that satisfy the one sd rule
# then its minimum value will be the value of d chosen by cross validation under ml

if (obj.mv[(dmax+1)] <= min(obj.mv) + obj.mv.sd[indice.min.mv])
{dselec.mv<-c(dselec.mv,0)}
#if indice.min.mv == (dmax+1), previous if condition will trivially be TRUE, because sd > 0, so dselec.mv will be equal
#to c((dmax+1),0) and finally, d.crossval.mv<-min(dselec.mv)==0


if (indice.min.mv > 1){
  for (dd in (1:(indice.min.mv-1))){
    if (obj.mv[dd] < min(obj.mv) + obj.mv.sd[indice.min.mv])
      dselec.mv<-c(dselec.mv,dd)}}

# if (indice.min.mv == 1) then nothing else has to be done, because comparison with d=0 has been done in the first step

d.crossval.mv<-min(dselec.mv)

#same for tau estimation
indice.min.rob<-which.min(funcion.FiB)  # is a number between 1 and 6
dselec.rob<-indice.min.rob

if (funcion.FiB[(dmax+1)] <= min(funcion.FiB) + funcion.FiB.sd[indice.min.rob])
{dselec.rob<-c(dselec.rob,0)}
#if indice.min.mv == (dmax+1), previous if condition will trivially be TRUE, because sd > 0, so dselec.mv will be equal
#to c((dmax+1),0) and finally, d.crossval.mv<-min(dselec.mv)==0


if (indice.min.rob > 1){
  for (dd in (1:(indice.min.rob-1))){
    if (funcion.FiB[dd] < min(funcion.FiB) + funcion.FiB.sd[indice.min.rob])
      dselec.rob<-c(dselec.rob,dd)}}

d.crossval.rob<-min(dselec.rob)


# now we plot
if (grafico=="TRUE"){
  #win.graph()
  graphics::par(mfrow=c(1,2))

  centro<-obj.mv
  desvio<-obj.mv.sd
  Cen =c(centro[(dmax+1)],centro[1:dmax])

  centror<-funcion.FiB
  desvior<-funcion.FiB.sd
  Cenr =c(centror[(dmax+1)],centror[1:dmax])


  #require(plotrix)
  plotrix::plotCI(0:dmax,  Cen, ui=Cen+c(desvio[(dmax+1)],desvio[1:dmax]), li=Cen-c(desvio[(dmax+1)],desvio[1:dmax]),
         xlab="d",
         ylab = "",main="MLE")
  graphics::lines(0:dmax, c(centro[(dmax+1)],centro[1:dmax]),col="blue")
  graphics::mtext(side=2,text=expression(paste(ln,"(",L,paste("(",hat(theta[d]),")",")"))), line = 2)
  if ( d.crossval.mv == 0 ){
    graphics::points(x = d.crossval.mv, y = obj.mv[(dmax+1)],col="red",pch=16)
    } else{graphics::points(x = d.crossval.mv, y = obj.mv[d.crossval.mv],col="red",pch=16) }

  plotrix::plotCI(0:dmax, Cenr , ui=Cenr+c(desvior[(dmax+1)],desvior[1:dmax]), li=Cenr -c(desvior[(dmax+1)],desvior[1:dmax]),
         xlab="d",ylab="",
        main="Tau")
  graphics::lines(0:dmax, c(centror[(dmax+1)],centror[1:dmax]),col="blue")
  graphics::mtext(side = 2, text = expression(paste(ln,"(",Phi[B],paste("(",hat(theta[d]),")",")"))), line = 2)
  if ( d.crossval.rob == 0 ){
    graphics::points(x=d.crossval.rob,y=funcion.FiB[(dmax+1)],col="red",pch=16)
  }  else{graphics::points(x=d.crossval.rob,y=funcion.FiB[d.crossval.rob],col="red",pch=16)}

  graphics::par(mfrow=c(1,1))

}

list(d.crossval.ml=d.crossval.mv,obj.ml=obj.mv[c((dmax+1),1:dmax)],obj.ml.sd=obj.mv.sd[c((dmax+1),1:dmax)],
     d.crossval.rob=d.crossval.rob,
     obj.rob=funcion.FiB[c((dmax+1),1:dmax)],obj.rob.sd=funcion.FiB.sd[c((dmax+1),1:dmax)])

}

##########################################################################################################
##########################################################################################################

initial_d0=function(X,aux,efficiency=0.85)
{ # INPUT
  # data matrix X
  # optionally the efficiency of the univariate estimate
  #
  # OUTPUT
  # beta0: an initial value of the coefficients
  # delta0: the S robust covariance as initial value of the covariance matrix
  #

  n=dim(X)[1]
  p=dim(X)[2]

  c1=aux$c1
  c2=aux$c2
  k1=aux$k1
  k2=aux$k2

  BETA=matrix(nrow=p,ncol=1)
  for (kkk in 1:p){

    BETA[kkk,]=robust::lmRob(X[,kkk]~1,control =robust::lmRob.control(efficiency = 0.85))$coefficients
  }
  aux2=X-matrix(rep(t(BETA),n),nrow=n,ncol=p,byrow=TRUE)

  A=rrcov::CovSest(aux2,method="bisquare")  #library(rrcov)
  delta0.sc=rrcov::getCov(A)

  dy=sqrt(diag((aux2)%*%solve(delta0.sc)%*%t(aux2)))
  s= MscaleNR(dy,c1,k1)$s

  dys=dy/s
  cte=mean(rhoc(dys,c2))*s^2/k2
  delta0=delta0.sc*cte
  list(beta0=BETA, delta0=delta0)

}

###############################################################################
tauestimate_d0=function(X,aux,inic)
{
  n=dim(X)[1]
  p=dim(X)[2]

  c1=aux$c1
  c2=aux$c2
  k1=aux$k1
  k2=aux$k2

  #####################################################
  delta0=inic$delta0
  beta0=inic$beta0

  ## Step 1
  aux2=X-matrix(rep(t(beta0),n),nrow=n,ncol=p,byrow=TRUE)
  dy=Re(sqrt(diag((aux2)%*%solve(delta0)%*%t(aux2))))
  s=MscaleNR(dy,c1,k1)$s
  dys=dy/s

  GG=NULL
  G0=phi1(delta0,c2,k2,s,dy,p)  # initial value of the function to be minimized
  GG[1]=G0

  w=weights2(dys,s,p,c1,c2)
  W=mean(w)

  W.matrix=matrix(rep(0,n^2),n,n)
  diag(W.matrix)=w
  unos=rep(1,n)
  mues=t(X)%*%W.matrix%*%unos/(n*W)

  aux2=X-matrix(mues,nrow=n,ncol=p,byrow=T)

  delta.sc=t(aux2)%*%W.matrix%*%aux2/n

  dy=Re(sqrt(diag((aux2)%*%solve(delta.sc)%*%t(aux2))))

  s=MscaleNR(dy,c1,k1)$s
  dys=dy/s
  cte=mean(rhoc(dys,c2))*s^2/k2
  deltaes=delta.sc*cte

  G1=phi1(deltaes,c2,k2,s,dy,p)

  GG[2]=G1

  k=2
  while ((abs(log(G0)-log(G1))>1e-5)& k<=20)
  {
    k=k+1

    ## Step 1
    aux2=X-matrix(rep(1,n),n,1)%*%t(mues)
    dy=Re(sqrt(diag((aux2)%*%solve(deltaes)%*%t(aux2))))
    s=MscaleNR(dy,c1,k1)$s
    dys=dy/s

    w=weights2(dys,s,p,c1,c2)
    W=mean(w)

    ##Step 2
    W.matrix=matrix(rep(0,n^2),n,n)
    diag(W.matrix)=w
    unos=rep(1,n)
    ##Step 3
    mues=(t(X)%*%W.matrix%*%unos)/(n*W)
    aux2=X-matrix(rep(1,n),n,1)%*%t(mues)

    delta.sc=t(aux2)%*%W.matrix%*%aux2/n


    ##Step 4
    dy=Re(sqrt(diag((aux2)%*%solve(delta.sc)%*%t(aux2))))

    s=MscaleNR(dy,c1,k1)$s
    dys=dy/s
    cte=mean(rhoc(dys,c2))*s^2/k2
    deltaes=delta.sc*cte
    dy=sqrt(diag((aux2)%*%solve(deltaes)%*%t(aux2)))
    G1=phi1(deltaes,c2,k2,s,dy,p)
    GG[k]=G1

  }

  alerta=(k>20)

  #-------
  list(mu=mues,delta=deltaes,alerta=alerta)
}
