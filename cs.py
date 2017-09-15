# timer:
from time import time
time0=time()
def pr(s,gn): 
      print gn,' :: ',time()-s
      return time(),gn+1
s,gn=pr(time0,0)

#import:
from numpy import pi, exp, log, sqrt,sin,cos
#from matplotlib.pyplot import *
from scipy.special import kn
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
#import numpy as np
#import matplotlib.pyplot as plt

#precomputations :
GaussLaguerre = laggauss(20) # laguerre-gauss sample points and weights
GaussLegendre = leggauss(40) # legendre-gauss sample points and weights
      # notice that if the numer of points is even then there is a zero angle or a cosine which equals to 1
      # That will lead to Zero Division Error in this function (sometimes)


s,gn=pr(s,gn)
def sigma_cs(x): # Mean compton scattering cross-section for electron rest frame i belive
      if x<.1:
            a,n,s=3./8 , 0. , 0.
            while(abs(a)*(n+2)**2>1e-11): # Tailor series sum of the formula below
                  s=s+a*(n+2+2/(n+1)+8/(n+2)-16/(n+3))
                  n=n+1
                  a=-2*x*a
            return s
      else: return 3*(2-(1/x+1-x/2)*log(1+2*x))/4/x/x + 3*(1+x)/4/(1+2*x)**2


def Compton_redistribution_m(x1,x2,mu,gamma): 
      """Compton redistribution matrix for monoenergetic electron gas
      x1 and x2 - photon energies in units of electron rest mass ( h \nu / m_e c^2 ) 
      mu - cosine of a scattering angle 
      gamma - energy of each electron in the gas in units of the electron rest mass
      This fuctions returns (R,RI,RQ,RU)
      which are scalar redistribution functions for isotropic monoenergetic gas
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectfully
      R44 or RV is also not equal to zero but we never need it 
      """

      #the next variables' namesare adopted from J. Poutanen & O. Vilhu 1993
      r = (1.+mu)/(1.-mu)
      a1 = sqrt( (gamma-x1)**2 + r )
      a2 = sqrt( (gamma+x2)**2 + r )
      v = a1*a2
      u = a2-a1  #(x1+x2)*(2.*gamma+x2-x1)/(a1+a2) # the formulas probably give always the same numbers 
      
      q = x1*x2*(1.-mu)
      Q = sqrt( x1*x1 + x2*x2 - 2.*x1*x2*mu )
      gammaStar = (x1-x2+Q*sqrt( 1. + 2./q ) )/2.
#      print r,v,u,q,Q,gammaStar

      if gamma < gammaStar : 
            print 'w'
            return  (0.,0.,0.,0.)
      else:
            Ra = u*( u*u - Q*Q )*( u*u + 5.*v )/2./q/q/v/v/v + u*Q*Q/q/q/v/v
            Rb = 2./Q + u/v*( 1. - 2./q )
            Rc = u/q/v*( ( u*u - Q*Q )/r/q - 2.)
            Rd = 2./Q + 2*(u-Q)/r/q*((u-Q)/r/q*(2.*Q+u) - 4.) + 2.*u/v/q 

            R = Ra + Rb
            RI = Ra + Rc
            RU = Rd + 2.*Rc
            RQ = RU + Ra
            #if gamma- gammaStar < 2 : print ';;;',gamma- gammaStar

            #print r,v,u,q,Q,gammaStar,Ra,Rb,Rc,R

      
            return (R,RI,RQ,RU)




def Maxwell_r(Theta,gamma):
      """The normalized relativistic Maxwellian distribution
      the density of particles in the dimensionless momentum volume (4 \pi z^2 dz) is nomalized to unity
      Theta is the dimensionless electron gas temperature (Theta = k * T_e / m_e c^2)
      gamma is electron energy in units of the electron rest mass
      The fuction returns the momentum dencity value ( f(\gamma) )
      """
      # kn is a modified Bessel function
      r = 1./4./pi/Theta*exp(-gamma/Theta)/kn(2,1./Theta)
      return r

def Compton_redistribution(x1,x2,mu,Theta): 
      """thermal Compton redistribution matrix (integrated with electron distribution function)
      x1 and x2 - photon energies in units of electron rest mass ( h \nu / m_e c^2 ) 
      mu - cosine of a scattering angle 
      This fuctions returns (R,RI,RQ,RU)
      which are scalar redistribution functions for Maxwellian relativistic gas
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectfully 
      R44 or RV is also not equal to zero but we never need it  
      """

      q = x1*x2*(1.-mu)
      Q = sqrt( x1*x1 + x2*x2 - 2.*x1*x2*mu )
      gammaStar = (x1-x2+Q*sqrt( 1. + 2./q ) )/2.
      C=3./8.*Maxwell_r(Theta,gammaStar)
      
     
      point, weight = GaussLaguerre 
      NL = range(len(point)) # number of laguerre-gauss sample points
      
      summands = map(lambda i : map(lambda x: C*x*weight[i],\
           Compton_redistribution_m(x1,x2,mu,point[i]+gammaStar)),NL)
      R = [sum([summands[j][i] for j in NL]) for i in range(4)] # summing to obtain the integrals
      
      # # the next 5 lines do the same thing as previous three but faster for some reason
      # # ( actually, only 2% faster, not very helpful) 
      # R=[0.,0.,0.,0.]
      # for i in NL:
      #       T=Compton_redistribution_m(x1,x2,mu,point[i]+gammaStar)
      #       for j in range(4):
      #             R[j]+=C*weight[i]*T[j]

      return tuple(R)
      

def Compton_redistribution_aa(x1,mu1,x2,mu2,Theta):
      """azimuth-avereged Compton redistribution matrix 
      for computing of electron scattering source function 
      x1 and x2 are photon energies in units of electron rest mass ( h \nu / m_e c^2 ) 
      mu1 and mu2 are cosines of angles between photon propagation directions and fixed direction
      This function returns R11 R12 R21 R22 matrix elements
      We need only 2x2 matrix in the upper left corner of the general matrix,
      becouse U and V components on the Stokes vector are zero in this case.
      """
      
      mur1 = sqrt( 1. - mu1*mu1 ) # sinuses of the angles 
      mur2 = sqrt( 1. - mu2*mu2 ) # what the R stands for? I dunno.
      
      point, weight = GaussLegendre
      NL = range(len(point)) # number of laguerre-gauss sample points
      
      az_c=cos(pi*point)  # list of azimuth cosines
      az_s=sin(pi*point)  # list of azimuth sinuses
      sc_c=mu1*mu2-mur1*mur2*az_c # list of scattering angles cosines
      cos2chi1 = map(lambda c : 2.*(mu2-mu1*c)**2/mur1/mur1/(1.-c*c) - 1. , sc_c)
      cos2chi2 = map(lambda c : 2.*(mu1-mu2*c)**2/mur2/mur2/(1.-c*c) - 1. , sc_c)
      sin2chi1 = map(lambda c, s : -(mu1-mu2*c)*mur2/mur1/(1.-c*c) , sc_c , az_s)
      sin2chi2 = map(lambda c, s :  (mu1-mu2*c)*mur1/mur2/(1.-c*c) , sc_c , az_s)
      # cos2chi12 = map(lambda c,s : ((mu2-mu1*c)**2-mur1*mur1*mur2*mur2*s*s)/mur1/mur1/(1.-c*c) , sc_c,az_s)
      # print max(map(lambda a,b:abs(a-b)/(a+b), cos2chi12,cos2chi1)) # comparing two formulas. It turns out they're the same
      
      R=[[0.,0.],[0.,0.]]
      for i in NL:
            (C,I,Q,U)=Compton_redistribution(x1,x2,sc_c[i],Theta)
            R[0][0]+=C*pi*weight[i]
            R[0][1]+=I*pi*cos2chi2[i]*weight[i]
            R[1][0]+=I*pi*cos2chi1[i]*weight[i]
            R[1][1]+=(Q*pi*cos2chi1[i]*cos2chi2[i]+U*pi*sin2chi1[i]*sin2chi2[i])*weight[i]
      
      
      return R




s,gn=pr(s,gn)
x1,x2=2.,3.
for NL in range(1,10):
      mu1,mu2=.1,.2
      print Compton_redistribution_aa(x1,mu1,x2,mu2,.1/NL)

s,gn=pr(s,gn)

print 'Total time : ', s-time0
# ps, garbage=  laggauss(20) 
# plot(range(20),log(ps))

# def pr(x):
#       print x,sn(x)

# ars=map(lambda a: 1e2*2.**(-a),(range(20)))
# print map(sn,ars)

# plot(ars,map(lambda a: (sn(1e4*2.**(-a))),ars))

# yscale('log')
# xscale('log')
# show()
