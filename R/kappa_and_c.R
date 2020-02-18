kappa_and_c=function(p,efficiency)
{##-------------
  ##iNPUT:
  ## p=dimension vector
  ## efi=estimator efficiency
  ##OUTPUT:
  ## k1: value of E(rho1) to obtain breakdown point 0.5
  ## k2: value of E(rho2) to obtain efficiency (ARE) equal to efi
  ## c1 and c2: values for the tuning parameters for the quasi optimal functions defined
          #Salibian-Barrera & Yohai (2016)
   ##require(rhoc_f)



efi=efficiency
tablasi=(p<=20)&(efi==0.85|efi==0.90|efi==0.95)
if (!tablasi)
  {message("It is not an efficient procedure, so it can take some seconds")
integrand=function(x,c1,p)
{rhoc(abs(sqrt(x)),c1)*stats::dchisq(x,p)}

ff=function(c1,p)
{maximoderho<-rhoc(1200,c1)
stats::integrate(integrand,0,Inf,c1=c1,p=p)$value-0.5*maximoderho}
#####################################################
## auxiliar compute to the integration limits
x1=.1
x2=.11
paso=.01
while ((ff(x1,p)*ff(x2,p)>0)&(x2<=5*p))
{ x1=x2
  x2=x2+paso
}
inicialc1=x1
finalc1=x2
################################################
c1=stats::uniroot(ff,c(inicialc1,finalc1),p=p)$root
k1<-stats::integrate(integrand,0,Inf,c1=c1,p=p)$value


k0<-1
inteC=function(v,c2,k0,p)
(2*rhoc(sqrt(v)/k0,c2)-psic(sqrt(v)/k0,c2)*sqrt(v)/k0)*stats::dchisq(v,p)

inteD=function(v,c1,k0,p)
psic(sqrt(v)/k0,c1)*sqrt(v)/k0*stats::dchisq(v,p)

psiast=function(v,c2,c1,k0,p)
{
C=stats::integrate(inteC,0,Inf,c2=c2,k0=k0,p=p)$value
D=stats::integrate(inteD,0,Inf,c1=c1,k0=k0,p=p)$value
C*psic(v,c1)+D*psic(v,c2)
}

psiastprime=function(v,c2,c1,k0,p)
{
C=stats::integrate(inteC,0,Inf,c2=c2,k0=k0,p=p)$value
D=stats::integrate(inteD,0,Inf,c1=c1,k0=k0,p=p)$value
C*psicprime(v,c1)+D*psicprime(v,c2)
}

inte3=function(v,c1,c2,k0,p)
{
  ((p-1)*psiast(sqrt(v)/k0,c2,c1,k0,p)/(sqrt(v)/k0)+psiastprime(sqrt(v)/k0,c2,c1,k0,p))*stats::dchisq(v,p)
}

aux1=function(v,c2,c1,k0,p)
{psiast(sqrt(abs(v))/k0,c2,c1,k0,p)^2*stats::dchisq(v,p)}

aux2=function(v,p)
{sqrt(v)^2*stats::dchisq(v,p)}

fc2=function(c2,c1,k0,p)
{
den=stats::integrate(inte3,0,Inf,c1=c1,c2=c2,k0=k0,p=p)$value
c0=p/den
efi*c0^2*k0^2*stats::integrate(aux1,0,Inf,c2=c2,c1=c1,k0=k0,p=p)$value-stats::integrate(aux2,0,Inf,p=p)$value
}

#####################################################
## auxiliar compute to the integration limits
x1=.1
x2=.1+paso
while ((fc2(x1,c1,k0,p)*fc2(x2,c1,k0,p)>0)&(x2<=5*p))
{ x1=x2
  x2=x2+paso
}
inicialc2=x1
finalc2=x2
###############################################################
c2=stats::uniroot(fc2,c(inicialc2,finalc2),c1=c1,k0=k0,p=p)$root #  3.263612 for p=1!!

inte4=function(v,c2,p)
rhoc(abs(sqrt(v)),c2)*stats::dchisq(v,p)

k2=stats::integrate(inte4,0,Inf,c2=c2,p=p)$value

list( c1=c1, k1=k1, c2=c2, k2=k2)
###--------------
}
else {
if(efi==.85) {cat(" 1.212013  2.7563434 0.4999986 0.17456287
  2.080258  3.0193607 0.4999999 0.28976213
                   2.691675  3.2881063 0.4999960 0.36641554
                   3.184632  3.5219602 0.4999984 0.42542091
                   3.607571  3.7251514 0.4999984 0.47444313
                   3.983408  0.9251285 0.4999991 0.99880413
                   4.324898  1.5247912 0.5000026 0.99273576
                   4.640000  1.9392661 0.5000048 0.98588694
                   4.934024  2.2794256 0.4999970 0.97953706
                   5.210777  2.5762449 0.4999997 0.97386266
                   5.472946  2.8450772 0.4999993 0.96871474
                   5.722647  3.0922767 0.4999994 0.96408328
                   5.961571  3.3226418 0.4999998 0.95987653
                   6.190947  3.5403625 0.4999998 0.95594503
                   6.411914  3.7478566 0.4999996 0.95220128
                   6.625273  3.9459730 0.4999994 0.94867638
                   6.831870  4.1409175 0.4999996 0.94484142
                   7.032237  4.3288463 0.5000011 0.94113611
                   7.226928  4.5126706 0.4999995 0.93730965
                   7.416421  4.6948352 0.4999995 0.93314070",file="aa")
bb=matrix(scan("aa",quiet = TRUE),nrow=20,byrow = TRUE)
ctes=bb[p,]}
  if (efi==.9){cat(" 1.212013 2.955710 0.4999986 0.1539176
2.080258 3.233177 0.4999999 0.2567185
                  2.691675 3.509092 0.4999960 0.3269221
                  3.184632 3.755970 0.4999984 0.3804467
                  3.607571 3.975865 0.4999984 0.4242064
                  3.983408 4.174260 0.4999991 0.4614332
                  4.324898 4.356999 0.5000026 0.4936464
                  4.640000 1.408707 0.5000048 0.9982838
                  4.934024 1.845705 0.4999970 0.9951378
                  5.210777 2.176719 0.4999997 0.9920144
                  5.472946 2.459318 0.4999993 0.9891192
                  5.722647 2.709519 0.4999994 0.9865737
                  5.961571 2.937468 0.4999998 0.9843326
                  6.190947 3.149394 0.4999998 0.9823181
                  6.411914 3.347599 0.4999996 0.9805440
                  6.625273 3.534795 0.4999994 0.9789591
                  6.831870 3.713259 0.4999996 0.9775057
                  7.032237 3.886090 0.5000011 0.9760617
                  7.226928 4.048204 0.4999995 0.9749425
                  7.416421 4.206679 0.4999995 0.9737876",file="aa")
bb=matrix(scan("aa",quiet = TRUE),nrow=20,byrow = TRUE)
ctes=bb[p,]}
    if (efi==.95){cat(" 1.212013 3.263646 0.4999986 0.1278930
2.080258 3.558998 0.4999999 0.2149413
                    2.691675 3.841485 0.4999960 0.2767617
                    3.184632 4.097121 0.4999984 0.3244408
                    3.607571 4.329569 0.4999984 0.3631703
                    3.983408 4.540652 0.4999991 0.3961612
                    4.324898 4.736932 0.5000026 0.4245809
                    4.640000 4.921284 0.5000048 0.4494459
                    4.934024 5.094868 0.4999970 0.4716094
                    5.210777 5.261084 0.4999997 0.4912821
                    5.472946 1.681907 0.4999993 0.9995862
                    5.722647 2.102544 0.4999994 0.9985270
                    5.961571 2.395138 0.4999998 0.9974579
                    6.190947 2.641974 0.4999998 0.9964377
                    6.411914 2.859370 0.4999996 0.9955312
                    6.625273 3.058209 0.4999994 0.9947152
                    6.831870 3.242833 0.4999996 0.9939916
                    7.032237 3.417637 0.5000011 0.9933164
                    7.226928 3.581752 0.4999995 0.9927453
                     7.416421 3.738933 0.4999995 0.9922213",file="aa")
bb=matrix(scan("aa",quiet = TRUE),nrow=20,byrow = TRUE)
ctes=bb[p,]}


list( c1=ctes[1], k1=ctes[3], c2=ctes[2], k2=ctes[4])
}
##
}
