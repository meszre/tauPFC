psic=function (x,c)    #derivative of rhoc
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
  v*G0*2*(x/c)/c+ w*(G2*2*(x/c)/c+G3*4*(x/c)^3/c+G4*6*(x/c)^5/c+ G5*8*(x/c)^7/c)+ 0*u
}