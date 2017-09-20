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
from numpy import linspace, empty, ones, zeros
from numpy import pi, exp, log, sqrt, sin, cos
from scipy.special import kn
#from matplotlib.pyplot import *
#import numpy as np
#import matplotlib.pyplot as plt
s,gn=pr(s,gn, 'importing')

#precomputations :
ScatterNum = 2 # total number of scatterings
NGamma=20 # number of Lorenz factor points (\gamma)
NAzimuth=40 # (*) numbers of azimuth angles (\phi)
NDirection = 10 # (*) number of propagation angle cosines (\mu) 
NEnergy = 11 # number of energy points (x)
NDepth = 23 # number of optical depth levels (\tau)
      # (*) Notice that if some numer of points is odd then there may be a zero angle or a cosine which equals to 1
      # That may lead to Zero Division Error in these functions, be careful and cautious 
tau_T= 1. # Thomson optical depth of thermalization 
x_l, x_u = -3. , 1. # lower and upper bound of the log_10 energy span 
logEnergySpan, logEnergyStep = linspace(x_l*log(10),x_u*log(10),num=NEnergy,retstep=True)

IntGamma = laggauss(NGamma) # sample points and weights for computing thermal matrix
IntAzimuth = leggauss(NAzimuth) # sample points and weights for computing azimuth-averaged matrix
IntDirection = leggauss(NDirection) #  sample points and weights for computing the source function of one given incoming photon energy
IntEnergy = (exp(logEnergySpan), exp(-logEnergySpan)*logEnergyStep,) # sample points and wights for computing the full source function 
IntDepth = linspace(0,tau_T,num=NDepth), tau_T/(NDepth-1)*ones(NDepth) # sample points and weights for computing the intensities

evere=.5109989e6 # electron volts / elecron rest energy 

def setTe(T):
      "use this function to set electron gas Maxwellian temperature"
      K = kn(2,1./T) # second modified Bessel function of reversed dimensionless temperature       
      return (T, K) # Theta is the dimensionless electron gas temperature (Theta = k * T_e / m_e c^2)

Theta, K2Y = setTe(5e4/evere) # it's about 0.1
T = 1e3/evere # photon black body temperature
s,gn=pr(s,gn,'precomps')


def Planck(x):
      """Planck function for Intensity of black body radiation 
      x is the energy of a photon in units of electron rest energy ( h \nu / m_e c^2 ) 
      The photon temperature is given by T also in units of electron rest mass
      Planck returns the intensity of  BB radiation
      """
      eps=1e-5
      C=1. # some dimension constant. I just didn't think what it should be like yet
      R=C*x*x # Rayleigh Jeans law
      I=R*x/(exp(x/T)-1.) if x/T < eps else R*T
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
      """Compton redistribution matrix for monoenergetic electron gas
      x1 and x2 - photon energies in units of electron rest energy ( h \nu / m_e c^2 ) 
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
      """thermal Compton redistribution matrix (integrated with electron distribution function)
      And the distribution is maxwellian (if it's not the function must be modified)
      x1 and x2 - photon energies in units of electron rest energy ( h \nu / m_e c^2 ) 
      mu - cosine of a scattering angle 
      This fuctions returns (R,RI,RQ,RU)
      which are scalar redistribution functions for Maxwellian relativistic gas
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectfully 
      R44 or RV is also not equal to zero but we never need it  
      """
      print mu
      
      q = x1*x2*(1.-mu)
      print q
      print sqrt( 1. + 2./q )

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
      """azimuth-avereged Compton redistribution matrix 
      for computing of electron scattering source function 
      x1 and x2 are photon energies in units of electron rest energy ( h \nu / m_e c^2 ) 
      mu1 and mu2 are cosines of angles between photon propagation directions and fixed direction
      This function returns R11 R12 R21 R22 matrix elements
      We need only 2x2 matrix in the upper left corner of the general matrix,
      becouse U and V components on the Stokes vector are zero in this case.
      """
      # this function gives not the same result as its old fortran ancestor for some reason
      # the difference is a factor depending of x1 and x2, but the relations between different elements are alright


      mur1 = sqrt( 1. - mu1*mu1 ) # sinuses of the angles 
      mur2 = sqrt( 1. - mu2*mu2 ) # what the R stands for? I dunno.
      
      point, point_weight = IntAzimuth
      NL = range(NAzimuth) # indeces list
      
      az_c=cos(pi*point)  # list of azimuth cosines
      az_s=sin(pi*point)  # list of azimuth sinuses
      sc_c=mu1*mu2-mur1*mur2*az_c # list of scattering angles cosines
      sc_ss = 1. - sc_c**2 # list of squared sinuses
      cos2chi1 = -1. + 2.*(mu2-mu1*sc_c)**2/mur1/mur1/sc_ss  # list[ cos( 2 \chi_1 ) ]
      cos2chi2 = -1. + 2.*(mu1-mu2*sc_c)**2/mur2/mur2/sc_ss  # list[ cos( 2 \chi_2 ) ]
      sin2chiP = -4.*(mu1-mu2*sc_c)*(mu2-mu1*sc_c)*(az_s/sc_ss)**2  # list[ sin( 2 \chi_1 )*sin( 2 \chi_2 ) ]
      
      # # This is actually not faster and not slower but contains more characters 
      # cos2chi1 = map(lambda c : 2.*(mu2-mu1*c)**2/mur1/mur1/(1.-c*c) - 1. , sc_c) # list[ cos( 2 \chi_1 ) ]
      # cos2chi2 = map(lambda c : 2.*(mu1-mu2*c)**2/mur2/mur2/(1.-c*c) - 1. , sc_c) # list[ cos( 2 \chi_2 ) ]
      # sin2chiP = map(lambda c, s : -4.*(mu1-mu2*c)*(mu2-mu1*c)*(s/(1.-c*c))**2 , sc_c , az_s) # list[ sin( 2 \chi_1 )*sin( 2 \chi_2 ) ]
      
      # sin2chi1 = map(lambda c, s : -(mu1-mu2*c)*mur2*s/mur1/(1.-c*c) , sc_c , az_s) # list[ sin( 2 \chi_1 ) ]
      # sin2chi2 = map(lambda c, s :  (mu1-mu2*c)*mur1*s/mur2/(1.-c*c) , sc_c , az_s) # list[ sin( 2 \chi_2 ) ]
      # cos2chi12 = map(lambda c,s : ((mu2-mu1*c)**2-mur1*mur1*mur2*mur2*s*s)/mur1/mur1/(1.-c*c) , sc_c,az_s)
      # print max(map(lambda a,b:abs(a-b)/(a+b), cos2chi12,cos2chi1)) # comparing two formulas. It turns out they're the same
      


      R=zeros( (2,2,) )
      for i in NL:
            (C,I,Q,U)=Compton_redistribution(x1,x2,sc_c[i])
            R[0][0]+=C*pi*point_weight[i]
            R[0][1]+=I*pi*cos2chi2[i]*point_weight[i]
            R[1][0]+=I*pi*cos2chi1[i]*point_weight[i]
            R[1][1]+=pi*(Q*cos2chi1[i]*cos2chi2[i]-U*sin2chiP[i])*point_weight[i]        
      
      
     # print x1,x2,mu1,mu2,R
      return R

s,gn=pr(s,gn,'funcs')


def CheckAngularSymmetry(x1,x2,mu1,mu2):
      print '344444444444444444444444444444444444444444444444444444444'
      eps=1e-10
      one = Compton_redistribution_aa(x1,x2,mu1,mu2)
      two = Compton_redistribution_aa(x1,x2,mu2,mu1)
      three = Compton_redistribution_aa(x1,x2,-mu1,-mu2)
      four = Compton_redistribution_aa(x1,x2,-mu2,-mu1) 
      v=True
      a,b,c = two[0][0]/one[0][0]-1,three[0][0]/one[0][0]-1,four[0][0]/two[0][0]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      a,b,c = two[1][0]/one[0][1]-1,three[1][0]/one[1][0]-1,four[0][1]/two[0][1]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      a,b,c = two[0][1]/one[1][0]-1,three[0][1]/one[0][1]-1,four[1][0]/two[1][0]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      a,b,c = two[1][1]/one[1][1]-1,three[1][1]/one[1][1]-1,four[1][1]/two[1][1]-1 # 0 0 0
      if abs(a)>eps or abs(b)>eps or abs(c)>eps : v=False
      return v

print CheckAngularSymmetry(0.01,0.1,-0.4,0.5)
s,gn=pr(s,gn,'ang-check')
exit()

for ComputingRedistributionMatrices in [1]:         
      RedistributionMatrix = empty( (NEnergy,NEnergy,NDirection,NDirection,2,2) )
      x,x_weight=IntEnergy
      mu,mu_weight=IntDirection
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

for InitializeStocksVectorsArrays in [1]:
      Source=empty((ScatterNum,NDepth,NEnergy,NDirection,2)) # source function                 
      Stokes=empty((ScatterNum,NDepth,NEnergy,NDirection,2)) # intensity Stokes vector
      Stokes_out=empty((ScatterNum+1,NEnergy,NDirection,2)) # outgoing Stokes vector of each scattering
      Stokes_in=empty((NDepth,NEnergy,NDirection,2)) # Stokes vector of the initial raiation (0th scattering) 
      Intensity=empty((NEnergy,NDirection,2)) # total intensity of all scattering orders from the slab suface 

      tau,tau_weight=IntDepth

      for e in range(NEnergy):
            for d in range(NDirection):
                  for t in range(NDepth):
                        Stokes_in[t][e][d][0]=Planck(x[e])*exp(-tau[t]/mu[d]) if mu[d]>0 else 0 
                        Stokes_in[t][e][d][1]=0
                  else:
                        Stokes_out[0][e][d][0]=Planck(x[e])*exp(-tau_T/mu[d]) if mu[d]>0 else 0
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
                                    w = x_weight[e1]*mu_weight[d1]
                                    r = RedistributionMatrix[e][e1][d][d1]
                                    I = Stokes[k-1][t][e1][d1] if k>0 \
                                          else Stokes_in[t][e1][d1]
                                    # print w, r, I
                                    S[0]+= w*( I[0]*r[0][0] + I[1]*r[0][1] )
                                    S[1]+= w*( I[0]*r[1][0] + I[1]*r[1][1] )
                        Source[k][t][e][d]=S*x[e]**2
                        # print 'S: ', S*x[e]**2
      
      for e in range(NEnergy): # I_k= integral S_k
            sigma=sigma_cs(x[e])
            for d in range(NDirection): 
                  m=mu[d]
                  for t in range(NDepth):
                        range_t = range(t) if m>0 else range (t+1,NDirection)
                        I=zeros(2)
                        for t1 in range_t:
                              w=tau_weight[t1]
                              S=Source[k][t1][e][d]
                              I+=w*S*exp(sigma*(tau[t1]-tau[t])/m)/m
                              # print I,w
                        Stokes[k][t][e][d]=I
                        # print 'I: ',I
                  else: 
                        Stokes_out[k+1][e][d]=Stokes[k][t][e][d]
      Intensity+=Stokes_out[k+1]       
      s,gn=pr(s,gn,'I'+str(1+k)) # do ScatterNum scattering iterations
print mu[NDirection/2:]
                                                

#print RedistributionMatrix




# ars=map(lambda a: 1e2*2.**(-a),(range(20)))
# print map(sn,ars)

# plot(ars,map(lambda a: (sn(1e4*2.**(-a))),ars))

# yscale('log')
# xscale('log')
# show()



# # frequency symmetry: CHECK [v]

one = Compton_redistribution_aa(1e-2,-.5,1e-1 ,.5)
two = Compton_redistribution_aa(1e-1,-.5,1e-2 ,.5) 
print exp(.09/Theta)*two[0][0]/one[0][0]-1 # 1e-15
print exp(.09/Theta)*two[1][0]/one[1][0]-1 # 1e-12
print exp(.09/Theta)*two[0][1]/one[0][1]-1 # 1e-12
print exp(.09/Theta)*two[1][1]/one[1][1]-1 # 1e-14
# s,gn=pr(s,gn,'freq-check')

# # angular symmetry: CHECK [x]


print 'Total time : ', s-time0
