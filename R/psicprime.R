psicprime=function (x,c)    #second derivative of rhoc
{
  G0=1.38
  G1=0.55
  G2=-2.69
  G3=10.76
  G4=-11.66
  G5=4.04
  u =(x>c)
  v=(x<2/3*c)
  w=(1-u)*(1-v)
  v*G0*2*(1/c^2)+ w*(G2*2*(1/c^2)+G3*4*3*(x/c)^2/c^2+G4*6*5*(x/c)^4/c^2+ G5*8*7*(x/c)^6/c^2)+ 0*u
}
