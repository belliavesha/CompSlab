# timer:
from time import time
time0=time()
def pr(s,gn,message = ''): 
      print gn,' :'+message+': ',time()-s
      return time(),gn+1
s,gn=pr(time0,0)

#import:
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from numpy import linspace, logspace, empty, zeros
from numpy import pi, exp, log, sqrt, sin, cos
from scipy.special import kn
from matplotlib.pyplot import *
#import numpy as np
#import matplotlib.pyplot as plt
s,gn=pr(s,gn, 'importing')

#physical constants:
evere=.5109989e6 # electron volts / elecron rest energy 

#precomputations :
ScatterNum = 10 # total number of scatterings
NGamma=10 # number of Lorenz factor points (\gamma)
NAzimuth=10 # (*) numbers of azimuth angles (\phi)
NDirection = 12 # (*) number of propagation angle cosines (\mu) 
NEnergy = 21 # number of energy points (x)
NDepth = 20 # number of optical depth levels (\tau)
      # (*) Notice that if some numer of points is odd then there may be a zero angle or a cosine which equals to 1
      # That may lead to Zero Division Error in these functions, be careful and cautious 
tau_T= .5 # Thomson optical depth of thermalization 
x_l, x_u = -6. , 0. # lower and upper bounds of the log_10 energy span

IntGamma = laggauss(NGamma) # sample points and weights for computing thermal matrix
IntAzimuth = leggauss(NAzimuth) # sample points and weights for computing azimuth-averaged matrix
IntDirection = leggauss(NDirection) #  sample points and weights for computing the source function of one given incoming photon energy
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for computing the full source function 
IntDepth = linspace(0,tau_T,num=NDepth,retstep=True)  # sample points and weights for computing the intensities

def setTe(T):
      " Use this function to set electron gas Maxwellian temperature"
      K = kn(2,1./T) # second modified Bessel function of reversed dimensionless temperature       
      return (T, K) # Theta is the dimensionless electron gas temperature (Theta = k * T_e / m_e c^2)

Theta, K2Y = setTe(0.1) # it's about 0.1
T = 1e1/evere # photon black body temperature
s,gn=pr(s,gn,'precomps')


def Planck(x):
      """   Planck function for Intensity of black body radiation 
      The only argument x is the energy of a photon in units of electron rest energy ( h \\nu / m_e c^2 ) 
      The photon temperature is given by T also in units of electron rest mass
      Planck returns the intensity of  BB radiation
      """
      eps=1e-5
      C=1.*evere # some dimension constant. I just didn't think what it should be like yet
      R=C*x*x # Rayleigh Jeans law
      I=R*x/(exp(x/T)-1.) if x/T > eps else R*T
      return I

def sigma_cs(x): # if Theta isn't small one will need to compute the mean cross-section  based on electron momentum distribution
      """ This function compute the Compton scattering cross-section in electron rest frame 
      x is the energy of the scattering photon in units of electron rest energy
      this function approaches the mean compton cross-section when electron gas temperature is small
      """ 
      if x<.1:
            a,n,s=3./8 , 0. , 0.
            while(abs(a)*(n+2)**2>1e-11): # Tailor series sum of the formula below
                  s=s+a*(n+2+2/(n+1)+8/(n+2)-16/(n+3))
                  n=n+1
                  a=-2*x*a
            return s
      else: return 3*(2-(1/x+1-x/2)*log(1+2*x))/4/x/x + 3*(1+x)/4/(1+2*x)**2

def Compton_redistribution_m(x1,x2,mu,gamma): 
      """   Compton redistribution matrix for monoenergetic electron gas
      The arguements are:
      x1 and x2 - photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
      mu - cosine of a scattering angle 
      gamma - energy of each electron in the gas in units of the electron rest mass
      
      This fuctions returns (R,RI,RQ,RU)
      which are scalar redistribution functions for isotropic monoenergetic gas
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectfully
      R44 or RV is also not equal to zero but we never need it 
      """

      #the next variables' names are adopted from J. Poutanen & O. Vilhu 1993
      r = (1.+mu)/(1.-mu)
      a1 = sqrt( (gamma-x1)**2 + r )
      a2 = sqrt( (gamma+x2)**2 + r )
      v = a1*a2
      u = a2-a1  #(x1+x2)*(2.*gamma+x2-x1)/(a1+a2) # the formulas probably give always the same numbers 
      
      q = x1*x2*(1.-mu)
      Q = sqrt( x1*x1 + x2*x2 - 2.*x1*x2*mu ) # sqrt( (x1-x2)**2 +2.*q ) # the formula probably gives the same also
      gammaStar = (x1-x2+Q*sqrt( 1. + 2./q ) )/2.

      if gamma < gammaStar : 
            print 'w' # I belive that in the case fucntion just won't be called
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
            
            #print r,v,u,q,Q,gammaStar,Ra,Rb,Rc,R
            #print x1,x2,mu,gamma,R,RI,RQ,RU
            
            return (R,RI,RQ,RU)
            #the return values of this functon seem to be all right

def Maxwell_r(gamma):
      """The normalized relativistic Maxwellian distribution
      the density of particles in the dimensionless momentum volume (4 \pi z^2 dz) is nomalized to unity
      Theta is the dimensionless electron gas temperature (Theta = k * T_e / m_e c^2)
      gamma is electron energy in units of the electron rest mass
      The fuction returns the momentum dencity value ( f(\gamma) )
      """
      r = .25/pi/Theta*exp(-gamma/Theta)/K2Y
      return r

def Compton_redistribution(x1,x2,mu): # if distribution is not Maxwellian the function must be modified.
      """    Thermal Compton redistribution matrix (integrated with electron distribution function)
      And the distribution is maxwellian (if it's not the function must be modified)
      The arguments are:
      x1 and x2 - photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
      mu - cosine of a scattering angle 
      
      This fuctions returns (R,RI,RQ,RU)
      which are scalar redistribution functions for Maxwellian relativistic gas
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectfully 
      R44 or RV is also not equal to zero but we never need it  
      """
      q = x1*x2*(1.-mu)
      Q = sqrt( x1*x1 + x2*x2 - 2.*x1*x2*mu )
      gammaStar = (x1-x2+Q*sqrt( 1. + 2./q ) )/2.
      C=3./8.*Theta*Maxwell_r(gammaStar)
            
     
      gamma, gamma_weight = IntGamma 
      NL = range(NGamma) # list of indeces
      
      summands = map(lambda i : map(lambda x: C*x*gamma_weight[i],\
           Compton_redistribution_m(x1,x2,mu,Theta*gamma[i]+gammaStar)),NL)
      R = [sum([summands[j][i] for j in NL]) for i in range(4)] # summing to obtain the integrals
      
      # # the next 5 lines do the same thing as previous three but faster for some reason
      # # ( actually, only 2% faster, not very helpful) 
      # R=[0.,0.,0.,0.]
      # for i in NL:
      #       T=Compton_redistribution_m(x1,x2,mu,Theta*gamma[i]+gammaStar)
      #       for j in range(4):
      #             R[j]+=C*gamma_weight[i]*T[j]

      return tuple(R)
      
def Compton_redistribution_aa(x1,x2,mu1,mu2):
      """   Azimuth-avereged Compton redistribution matrix 
      for computing of electron scattering source function 
      The arguements are:
      x1 and x2 are photon energies in units of electron rest energy ( h \\nu / m_e c^2 ) 
      mu1 and mu2 are cosines of angles between photon propagation directions and fixed direction
      
      This function returns R11 R12 R21 R22 matrix elements
      We need only 2x2 matrix in the upper left corner of the general matrix,
      becouse U and V components on the Stokes vector are zero in this case.
      """
      # this function gives not the same result as its old fortran ancestor for some reason
      # the difference is a factor depending of x1 and x2, but the relations between different elements are alright


      eta1 = 1. - mu1*mu1  # squared sinuses of the angles 
      eta2 = 1. - mu2*mu2  
      
      point, point_weight = IntAzimuth
      NL = range(NAzimuth) # indeces list
      
      az_c=cos(pi*point)  # list of azimuth cosines
      az_s=sin(pi*point)**2  # list of azimuth square sinuses
      sc_c=mu1*mu2-sqrt(eta1*eta2)*az_c # list of scattering angles' cosines
      sc_s=1. - sc_c**2 # list of scattering angles' squared sinuses
      cos2chi1 = 2.*(mu1*sc_c-mu2)*(mu1*sc_c-mu2)/eta1/sc_s-1.  # list[ cos( 2 \chi_1 ) ]
      cos2chi2 = 2.*(mu1-mu2*sc_c)*(mu1-mu2*sc_c)/eta2/sc_s-1.  # list[ cos( 2 \chi_2 ) ]
      sin2chiP = 4.*(mu1-mu2*sc_c)*(mu1*sc_c-mu2)*az_s/sc_s**2  # list[ sin( 2 \chi_1 )*sin( 2 \chi_2 ) ]

      R=zeros( (2,2,) )
      for i in NL:
            (C,I,Q,U)=Compton_redistribution(x1,x2,sc_c[i])
            R[0][0]+=C*pi*point_weight[i]
            R[0][1]+=I*pi*cos2chi2[i]*point_weight[i]
            R[1][0]+=I*pi*cos2chi1[i]*point_weight[i]
            R[1][1]+=pi*(Q*cos2chi1[i]*cos2chi2[i]+U*sin2chiP[i])*point_weight[i]         
      
      
     # print x1,x2,mu1,mu2,R
      return R*x1*x1/x2

s,gn=pr(s,gn,'funcs')

if True: #draw Planck
      x=IntEnergy[0]
      plot(x,map(lambda e: log(e*Planck(e)),x),'g')
      
if False: # compare the numbers to ones obtained from old fortran code
      from compare import Factor
      Factor(Compton_redistribution_aa)
      
if False: # checking symmetries
      from check import *
      print CheckAngularSymmetry(Compton_redistribution_aa,0.01,0.1,-0.4,0.5)
      s,gn=pr(s,gn,'ang-check')
      print CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.1,0.5,-0.5,Theta)
      s,gn=pr(s,gn,'freq-check')

if True: # Computing Compton scattering cross section for all energies
      x,x_weight=IntEnergy
      sigma=map(sigma_cs,x)   

if True : # Computing redistribution matrices for all energies and angles 
      mu,mu_weight=IntDirection
      RedistributionMatrix = empty( (NEnergy,NEnergy,NDirection,NDirection,2,2) )
      for e in range(NEnergy):
            for e1 in range(e,NEnergy):
                  for d in range(NDirection/2): 
                        for d1 in range(d,NDirection/2):
                              r=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[d1])

                              t=d1>d
                              f=e1>e
                 
                              md=NDirection-d-1 
                              md1=NDirection-d1-1

                              RedistributionMatrix[e][e1][d][d1]=r
                              RedistributionMatrix[e][e1][md][md1]=r
                              if t: #angular symmetry
                                    rt=r
                                    rt[0][1],rt[1][0]=rt[1][0],rt[0][1]
                                    RedistributionMatrix[e][e1][d1][d]=rt
                                    RedistributionMatrix[e][e1][md1][md]=rt
                              if f: #frequency symmetry
                                    m=exp((x[e]-x[e1])/Theta)
                                    rf=r*m
                                    RedistributionMatrix[e1][e][d][d1]=rf
                                    RedistributionMatrix[e1][e][md][md1]=rf
                              if t and f: # both
                                    rtf=rt*m
                                    RedistributionMatrix[e1][e][d1][d]=rtf
                                    RedistributionMatrix[e1][e][md1][md]=rtf
      s,gn=pr(s,gn,'RMtable')

if True : # Initializing Stokes vectors arrays, computing zeroth scattering 
      tau,tau_weight=IntDepth
      Source=empty((ScatterNum,NDepth,NEnergy,NDirection,2)) # source function                 
      Stokes=empty((ScatterNum,NDepth,NEnergy,NDirection,2)) # intensity Stokes vector
      Stokes_out=empty((ScatterNum+1,NEnergy,NDirection,2)) # outgoing Stokes vector of each scattering
      Stokes_in=empty((NDepth,NEnergy,NDirection,2)) # Stokes vector of the initial raiation (0th scattering) 
      Intensity=empty((NEnergy,NDirection,2)) # total intensity of all scattering orders from the slab suface 
      for e in range(NEnergy):
            for d in range(NDirection):
                  for t in range(NDepth):
                        Stokes_in[t][e][d][0]=Planck(x[e])*exp(-tau[t]*sigma[e]/mu[d]) if mu[d]>0 else 0 
                        Stokes_in[t][e][d][1]=0
                  else:
                        Stokes_out[0][e][d][0]=Planck(x[e])*exp(-tau_T*sigma[e]/mu[d]) if mu[d]>0 else 0
                        Stokes_out[0][e][d][1]=0
      Intensity=Stokes_out[0]
      s,gn=pr(s,gn,'I0')

for k in range(ScatterNum): # do ScatterNum scattering iterations
      for t in range(NDepth): # S_k= R T_{k-1}
            for e in range(NEnergy):
                  for d in range(NDirection):
                        S=zeros(2)
                        for e1 in range(NEnergy):
                              for d1 in range(NDirection):
                                    w = mu_weight[d1]*x_weight # total weight
                                    r = RedistributionMatrix[e][e1][d][d1]
                                    I = Stokes[k-1][t][e1][d1] if k>0 \
                                          else Stokes_in[t][e1][d1]
                                    # print w, r, I
                                    S[0]+= w*( I[0]*r[0][0] + I[1]*r[0][1] )
                                    S[1]+= w*( I[0]*r[1][0] + I[1]*r[1][1] )
                        Source[k][t][e][d]=S
                        # print 'S: ', S*x[e]**2
      
      
      for t in range(NDepth):
            for e in range(NEnergy): # I_k= integral S_k
                  for d in range(NDirection): 
                        # m=mu[d]
                        I=zeros(2)
                        for t1 in range(t) if mu[d]>0 else range (t+1,NDirection):
                              S=Source[k][t1][e][d]
                              I+=tau_weight*S*exp(sigma[e]*(tau[t1]-tau[t])/mu[d])
                              # print I,w
                        Stokes[k][t][e][d]=I/mu[d]
                        # print 'I: ',I
      else:
            Stokes_out[k+1]=Stokes[k][t]
            Intensity+=Stokes[k][t]       
      s,gn=pr(s,gn,'I'+str(1+k)) # do ScatterNum scattering iterations

out = open('res','w')
d=NDirection*3/5
out.write(str(mu[d])+str(list(x))+'\n')
#range(NDirection/2,NDirection):
print d
for k in range(ScatterNum+1):
      xFx=[log(Stokes_out[k][e][d][0]*x[e]) for e in range(NEnergy)]
      plot(x,xFx,'r')
      out.write(str(k)+' : '+str(mu[d])+str(xFx)+'\n')

xFx=[log(Intensity[e][d][0]*x[e]) for e in range(NEnergy)]
plot(x,xFx)

# diff = (lambda a,b : (log(Intensity[a][d][0])-log(Intensity[b][d][0]))/(log(x[a])-log(x[b])))
# print diff(10,12),diff(12,14) ,diff(14,16),diff(16,18)  
# yscale('log')
xscale('log')
out.write('t : '+str(mu[d])+str(xFx)+'\n')
show()                                           

print 'Total time : ', s-time0

