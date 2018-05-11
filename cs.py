
# coding: utf-8

# In[17]:



# timer:
from time import time
time0=time()
def pr(s,gn,message = ''): 
      print(gn,' :'+message+': ',time()-s)
      return time(),gn+1
s,gn=pr(time0,0)

Spectrum={
      1:'Compton',
      2:'Burst',
      3:'Thomson',
      0:'FromFile'
}[0]

oblateness='AlGendy'

AtmName='res/C/C3' # the prefix for all result files related to the set of parameters
PulsName=AtmName+'P12'
computePulse= True
plotAtm=True
plotPulse=True


# In[18]:



#import:
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.special import kn
from matplotlib.pyplot import *
from bisect import bisect

#import numpy as np
#import matplotlib.pyplot as plt
s,gn=pr(s,gn, 'importing')


# In[19]:




colors=['xkcd:brownish red',
        'xkcd:red',
        'xkcd:orange',
        'xkcd:dark yellow',
        'xkcd:dark yellow green',
        'xkcd:deep green',
        'xkcd:dark cyan',
        'xkcd:blue',
        'xkcd:purple'   
] # 'rygcbm' # Rainbow 
NColors=len(colors)

#physical constants:
evere=.5109989e6 # electron volts in elecron rest energy 
G=13275412528e1 # G*M_sol in km^3/s^2 
c=299792458e-3 # speed of light in km/s

# Atmosphere parameters: 
tau_T= .6 # Thomson optical depth of thermalization 
x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
Theta = 0.1  # dimensionless electron gas temperature (Theta = k T_e / m_e c^2) # it's about 0.1 
T = 0.002 # 10/evere #  dimensionless photon black body temperature T = k T_bb / m_e c^2

#precomputations :
ScatterNum = 23 # total number of scatterings
NGamma= 7# number of Lorenz factor points (\gamma)
NAzimuth= 10 # 12 # numbers of azimuth angles (\phi) [0,pi]
NEnergy = 281 # 50# 101 # number of energy points (x)
NDepth = 100# 101  # number of optical depth levels (\tau)
NMu = 22 # 20# 15 # number of propagation zenith angle cosines (\mu) [0,1]
NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]

IntGamma = laggauss(NGamma) # sample points and weights for computing thermal matrix
IntAzimuth = leggauss(NAzimuth*2) # sample points and weights for computing azimuth-averaged matrix
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function
IntDepth = linspace(0,tau_T,num=NDepth,retstep=True)  # sample points and weights for integrations over the optical depth computing intencity 
IntZenith = leggauss(NZenith) #  sample points and weights for integrations over zenith angle in positive and negative directions together

# IntZenith = linspace(1e-6-1,1-1e-6,num=NZenith), array([2.]+[2.]*(NZenith-2)+[2.])/(NZenith)
# IntZenith = linspace(1/NZenith - 1,1 - 1/NZenith,num=NZenith), array([2.]*(NZenith))/(NZenith)
# mu,mu_weight=IntZenith
# mu=array( ([1-1e-7]+list(mu)+[1-1e-7]) )
# mu_weight=array([1e-5]+list(mu_weight)+[1e-5])
# IntZenith=(mu, mu_weight)
# NZenith=NZenith+2
# NMu=NMu+1


K2Y = kn(2,1./Theta) # second modified Bessel function of reversed dimensionless temperature       

mu,mu_weight=IntZenith
x,x_weight=IntEnergy
tau,tau_weight=IntDepth

outParams = open(AtmName+'.dat','w')
outParams.write(Spectrum+'\n')
outParams.write(str(tau_T)+' = total depth tau_T\n')
outParams.write(str(x_l)+' '+str(x_u)+' = log_10 energy span from to\n')
outParams.write(str(Theta)+' = dimensionless electron gas temperature in m_ec^2\n')
outParams.write(str(T)+' = dimensionless initial bb photons temperature in m_ec^2\n')
outParams.write(str(ScatterNum)+' = Number of scatterings computed\n')
outParams.write(str(NMu)+' = NMu, NZenith=2 NMu, number of points in mu grid\n')
outParams.write(str(NEnergy)+' = NEnergy, number of points in x grid\n')
outParams.write(str(NGamma)+' '+str(NAzimuth)+' '+str(NDepth)+' = NGamma NAzimuth NDepth parameters\n')



s,gn=pr(s,gn,'precomps')


# In[20]:




def Planck(x):
      """   Planck function for Intensity of black body radiation 
      The only argument x is the energy of a photon in units of electron rest energy ( h \\nu / m_e c^2 ) 
      The photon temperature is given by T also in units of electron rest mass
      Planck returns the intensity of  BB radiation
      """
      e=x/T
      C=1.  # some dimension constant.
      R=C*e*e # Rayleigh Jeans law'
      I=.0 if e>5e2 else R*e/(exp(e)-1.) if e > 1e-5 else R
      return I

def Delta(x):
      C=2e4
      I=C*exp(-1e2*(x-T)**2/T/T)
      return I      

def sigma_cs(x): # not averaged on electron distribution 
      """ This function compute the Compton scattering cross-section in electron rest frame 
      x is the energy of the scattering photon in units of electron rest energy
      this function approaches the mean compton cross-section when electron gas temperature is small
      """ 
      if x<.1:
            a,n,s=3./8.,0.,0.
            while(abs(a)*(n+2)**2>1e-11): # Taylor series sum of the formula below
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
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectively 
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

      # print(  gamma-gammaStar, gammaStar,gamma,Q,u,(u-Q)/(Q+u))

      if gamma < gammaStar : 
            print ('w') # I belive that in the case fucntion just won't be called
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
            
            #print(r,v,u,q,Q,gammaStar,Ra,Rb,Rc,R)
            #print(x1,x2,mu,gamma,R,RI,RQ,RU)
            
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
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectively 
      R44 or RV is also not equal to zero but we never need it  
      """
      q = x1*x2*(1.-mu)
      Q = sqrt( x1*x1 + x2*x2 - 2.*x1*x2*mu )
      gammaStar = (x1-x2+Q*sqrt( 1. + 2./q ) )/2. # lower bound of integration 
      C=3./8.*Theta*Maxwell_r(gammaStar)

      gamma, gamma_weight = IntGamma 
      
      R=[0.,0.,0.,0.]
      for i in range(NGamma):
            T=Compton_redistribution_m(x1,x2,mu,Theta*gamma[i]+gammaStar)
            for j in range(4):
                  R[j]+=C*gamma_weight[i]*T[j]

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
      
      phi, phi_weight = IntAzimuth
      
      az_c=cos(pi*phi)  # array of azimuth cosines
      az_s=sin(pi*phi)**2  # array of azimuth square sinuses
      sc_c=mu1*mu2-sqrt(eta1*eta2)*az_c # array of scattering angles' cosines
      sc_s=1. - sc_c**2 # array of scattering angles' squared sinuses
      cos2chi1 = 2.*(mu1*sc_c-mu2)*(mu1*sc_c-mu2)/eta1/sc_s-1.  # array[ cos( 2 \chi_1 ) ]
      cos2chi2 = 2.*(mu1-mu2*sc_c)*(mu1-mu2*sc_c)/eta2/sc_s-1.  # array[ cos( 2 \chi_2 ) ]
      sin2chiP = 4.*(mu1-mu2*sc_c)*(mu1*sc_c-mu2)*az_s/sc_s**2  # array[ sin( 2 \chi_1 )*sin( 2 \chi_2 ) ]

      R=zeros( (2,2,),  )
      for i in range(NAzimuth*2):
            (C,I,Q,U)=Compton_redistribution(x1,x2,sc_c[i])
            R[0,0]+=C*pi*phi_weight[i]
            R[0,1]+=I*pi*cos2chi2[i]*phi_weight[i]
            R[1,0]+=I*pi*cos2chi1[i]*phi_weight[i]
            R[1,1]+=pi*(Q*cos2chi1[i]*cos2chi2[i]+U*sin2chiP[i])*phi_weight[i]         
      
      # print(x1,x2,mu1,mu2,R)
      return R*x1*x1/x2



# In[21]:



if Spectrum=='Compton' : # checking symmetries #unnecessary 
      from check import *
      print('Angular symmetry: ',CheckAngularSymmetry(Compton_redistribution_aa,0.02,0.02,-0.4,0.5))
      s,gn=pr(s,gn,'ang-check')
      print('Energy symmetry: ',CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.02,0.5,0.5,Theta))
      s,gn=pr(s,gn,'freq-check')
      # exit()


# In[22]:




if Spectrum=='Compton' : # Computing redistribution matrices for all energies and angles 
      # import FRM
      # print(FRM.__doc__)
      # print(FRM.fill.__doc__)
      # def cr(x1,x2,m1,m2):
      #       R=Compton_redistribution_aa(x1,x2,m1,m2)
      #       return 1.0,1.0,2.0,1.0
      #       # return R[0][0],R[0][1],R[1][0],R[1][1]
      # # FRM.cr=Compton_redistribution_aa
      # RedistributionMatrix , sigma = FRM.fill( mu_weight,x_weight,x,mu,Theta,cr)
      # print('7')
      # # print(sigma)
      # # exit()
      sigma=zeros(NEnergy)
      RedistributionMatrix = ones( (NEnergy,NEnergy,NZenith,NZenith,2,2) )
      percent=0.0
      for e in range(NEnergy): # x [-\infty,\infty]
            for e1 in range(e,NEnergy): # x1 [x,\infty]
                  percent+=200/NEnergy/(NEnergy+1)
                  npc=int(percent*0.6)
                  print('||'+'%'*npc+' '*(59-npc)+'|| {:5.3f}%'.format(percent))
                  for d in range(NMu): # mu [-1,0]
                        for d1 in range(d,NMu): # mu1 [-1,mu]
                              md=NZenith-d-1 # -mu
                              md1=NZenith-d1-1 # -mu1
                              w=mu_weight[d1]*x_weight*mu_weight[d]
                              t=d1>d
                              f=e1>e

                              r=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[d1])
                              rm=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[md1])
                              sigma[e1]+=(r[0,0]+rm[0,0])*w
                              RedistributionMatrix[e,e1,d,d1]=r
                              RedistributionMatrix[e,e1,md,md1]=r
                              RedistributionMatrix[e,e1,d,md1]=rm
                              RedistributionMatrix[e,e1,md,d1]=rm
                              if f: # frequency symmetry
                                    m=exp((x[e]-x[e1])/Theta)*x[e1]**3/x[e]**3
                                    rf=r*m  # when Maxwellian
                                    rmf=rm*m  # or Wein distributions
                                    sigma[e]+=(rf[0,0]+rmf[0,0])*w
                                    RedistributionMatrix[e1,e,d,d1]=rf
                                    RedistributionMatrix[e1,e,md,md1]=rf
                                    RedistributionMatrix[e1,e,d,md1]=rmf
                                    RedistributionMatrix[e1,e,md,d1]=rmf
                              if t: # angular symmetry
                                    r[0,1],r[1,0]=r[1,0],r[0,1]
                                    rm[0,1],rm[1,0]=rm[1,0],rm[0,1]
                                    sigma[e1]+=(r[0,0]+rm[0,0])*w
                                    RedistributionMatrix[e,e1,d1,d]=r
                                    RedistributionMatrix[e,e1,md1,md]=r
                                    RedistributionMatrix[e,e1,md1,d]=rm
                                    RedistributionMatrix[e,e1,d1,md]=rm
                                    if f: # both symmeties 
                                          rf=r*m
                                          rmf=rm*m
                                          sigma[e]+=(rf[0,0]+rmf[0,0])*w
                                          RedistributionMatrix[e1,e,d1,d]=rf
                                          RedistributionMatrix[e1,e,md1,md]=rf
                                          RedistributionMatrix[e1,e,md1,d]=rmf
                                          RedistributionMatrix[e1,e,d1,md]=rmf
      s,gn=pr(s,gn,'RMtable')





# In[26]:



if Spectrum=='Compton' : #check cross section  
      # cro = open('cros.dat','r')
      figure(1,figsize=(10,11))
      
      xscale('log')
      # fx=[0.0]
      # cs=[1.0]
      # sgm=[1.0]
      # de=lambda X : float (X[:X.find('D')]+'e'+X[X.find('D')+1:])
      # for n in range(85):
      #       line = cro.readline().split()
      #       fx.append(10**de(line[3]))
      #       cs.append(de(line[5]))
      #       sgm.append(sigma_cs(fx[-1]))
      #       # print x, cs
      # cs.append(0)
      # fx.append(x[-1])
      # sgm.append(sigma_cs(x[-1]))
      # fcs=interp1d(fx,cs)
      # fcsx=fcs(x)
      # plot(fx,cs)
      sigam=array(list(map(sigma_cs,x)))
      plot(x,sigam,'b')
      plot(x,sigma,'r')
      # plot(x,fcsx,'k')
      # plot(fx,sgm)
      # plot(x,sigma2-fcsx)
      print(sigma)
      print(sigam)
      print(sigam-sigma)
      print(x)
#       savefig('compsigma.png')
      # show()
      clf()
      s,gn=pr(s,gn,'sigmaplot') 


# In[27]:





if Spectrum=='Compton' : # Initializing Stokes vectors arrays, computiong scatterings 
      Iin=Planck # Delta # initial photon distribution 
      Source=zeros((ScatterNum,NDepth,NEnergy,NZenith,2)) # source function                 
      Stokes=zeros((ScatterNum,NDepth,NEnergy,NZenith,2)) # intensity Stokes vector
      Stokes_out=zeros((ScatterNum+1,NEnergy,NZenith,2)) # outgoing Stokes vector of each scattering
      Stokes_in=zeros((NDepth,NEnergy,NZenith,2)) # Stokes vector of the initial raiation (0th scattering) 
      Intensity=zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
      for e in range(NEnergy):
            for d in range(NZenith):
                  for t in range(NDepth):
                        Stokes_in[t,e,d,0]=Iin(x[e])*exp(-tau[t]*sigma[e]/mu[d]) if mu[d]>0 else 0 
                        Stokes_in[t,e,d,1]=0
                  else:
                        Stokes_out[0,e,d,0]=Iin(x[e])*exp(-tau_T*sigma[e]/mu[d]) if mu[d]>0 else 0
                        Stokes_out[0,e,d,1]=0
      s,gn=pr(s,gn,'I0')


      # if mod:
      try: # Fortran
            import FIF
            # print(FIF.fill.__doc__)
            Stokes=FIF.fill(ScatterNum,Stokes_in,RedistributionMatrix,x_weight,sigma,mu,mu_weight,tau,tau_weight)
            s,gn=pr(s,gn,'I')

      # else:
      except: # python
            for k in range(ScatterNum): # do ScatterNum scattering iterations
                  for t in range(NDepth): # S_k= R I_{k-1}
                        for e in range(NEnergy):
                              for d in range(NZenith):
                                    S=zeros(2)  
                                    for e1 in range(NEnergy):
                                          for d1 in range(NZenith):
                                                w = mu_weight[d1]*x_weight # total weight
                                                r = RedistributionMatrix[e,e1,d,d1]  #  
                                                I = Stokes[k-1,t,e1,d1] if k>0 else Stokes_in[t,e1,d1]
                                                S[0]+= w*( I[0]*r[0,0] + I[1]*r[0,1] ) # 
                                                S[1]+= w*( I[0]*r[1,0] + I[1]*r[1,1] ) #
                                    Source[k,t,e,d]+=S #     
                  
                  for t in range(NDepth):# I_k= integral S_k
                        for e in range(NEnergy): 
                              for d in range(NZenith): 
                                    I=tau_weight/2*Source[k,t,e,d]
                                    if mu[d]>0:
                                          for t1 in range(t) : #
                                                S=Source[k,t1,e,d] #
                                                I+=tau_weight*S*exp(sigma[e]*(tau[t1]-tau[t])/mu[d])
                                          S=Source[k,0,e,d] #
                                          I-=tau_weight*S*exp(sigma[e]*(-tau[t])/mu[d])
                                    else:
                                          for t1 in range (t+1,NDepth):
                                                S=Source[k,t1,e,d] #
                                                I+=tau_weight*S*exp(sigma[e]*(tau[t1]-tau[t])/mu[d])
                                          S=Source[k,NDepth-1,e,d] #
                                          I-=tau_weight*S*exp(sigma[e]*(tau_T-tau[t])/mu[d])
                                    Stokes[k,t,e,d]+=I/abs(mu[d]) #abs
                  s,gn=pr(s,gn,'I'+str(1+k))
                 
      Intensity += Stokes_out[0]
      for k in range(ScatterNum):
            Stokes_out[k+1]+=Stokes[k,-1]
            Intensity += Stokes[k,-1]


# In[ ]:





if Spectrum=='Thomson':
      Iin=Planck # Delta # initial photon distribution                  
      Stokes=zeros((ScatterNum,NDepth,NEnergy,NZenith,2)) # intensity Stokes vector
      Stokes_out=zeros((ScatterNum+1,NEnergy,NZenith,2)) # outgoing Stokes vector of each scattering
      Intensity=zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
      I_l = zeros((NDepth,NZenith))
      I_r = zeros((NDepth,NZenith))
      A=zeros(NDepth)
      B=zeros(NDepth)
      C=zeros(NDepth)
      x_factor=ones(NEnergy)
      x_shifted=x.copy()
      print('Thart')


      for d in range(NZenith):
            for t in range(NDepth):
                  I_l[t,d]=exp(-tau[t]/mu[d]) if mu[d]>0 else 0 
                  I_r[t,d]=I_l[t,d]
                  for e in range(NEnergy):
                        Stokes_out[0,e,d,0]=Iin(x[e])*(I_l[t,d] + I_r[t,d])
                        Stokes_out[0,e,d,1]=0
      for k in range(ScatterNum): # do ScatterNum scattering iterations
            for e in range(NEnergy):
                # print(k,e,.25 + 2*Theta + 4*Theta**2 - x_shifted[e])
                D=max(0,.25 + 2*Theta + 4*Theta**2 - x_shifted[e])
                #D=1 if x[e]> 4*Theta else 1+4*Theta
                x_shifted[e]=0.5+2*Theta-sqrt(D)
                x_factor[e] = 0 if x[e]> 4*Theta else x[e]/x_shifted[e] #cut off the high energies
            
            A[:]=0.
            B[:]=0.
            C[:]=0.
            for t in range(NDepth): #
                  for d in range(NMu):
                        md=NZenith-d-1
                        A[t]+=(1. - mu[d]**2)*(I_l[t,d] + I_l[t,md])*mu_weight[d]*3./4./3
                        B[t]+=(mu[md]**2)*(I_l[t,d] + I_r[t,md])*mu_weight[d]*3./8./3
                        C[t]+=(I_r[t,md] + I_r[t,d])*mu_weight[d]*3./8./3
            I_l[:]=0.
            I_r[:]=0.
            for t in range(NDepth): #
                  for d in range(NZenith):
                        w=tau_weight/abs(mu[d])
                        I_l[t,d]+=(A[t]*(1 - mu[d]**2) + (B[t] + C[t])*mu[d]**2)*w/2
                        I_r[t,d]+=(B[t] + C[t])*w/2
                        if mu[d]>0:
                              trange=range(t)
                              t1=0
                        else:
                              trange=range(t+1,NDepth)
                              t1=-1
                        I_l[t,d]-=exp((tau[t1]-tau[t])/mu[d])*(A[t1]*(1. - mu[d]**2) + (B[t1] + C[t1])*mu[d]**2)*w/2
                        I_r[t,d]-=exp((tau[t1]-tau[t])/mu[d])*(B[t1] + C[t1])*w/2
                        for t1 in trange:
                              I_l[t,d]+=exp((tau[t1] - tau[t])/mu[d])*(A[t1]*(1. - mu[d]**2) + (B[t1] + C[t1])*mu[d]**2)*w
                              I_r[t,d]+=exp((tau[t1] - tau[t])/mu[d])*(B[t1]+C[t1])*w

                        
                        for e in range(NEnergy):
                              Stokes[k,t,e,d,0]=Iin(x_shifted[e])*(I_l[t,d] + I_r[t,d])*x_factor[e]
                              Stokes[k,t,e,d,1]=Iin(x_shifted[e])*(I_l[t,d] - I_r[t,d])*x_factor[e]
            s,gn=pr(s,gn,'I'+str(k+1))


      Intensity += Stokes_out[0]
      for k in range(ScatterNum):
            Stokes_out[k+1]+=Stokes[k,-1]
            Intensity += Stokes[k,-1]
            

      s,gn=pr(s,gn,'I0')


# In[ ]:



if Spectrum=='Burst' : # Initializing Stokes vectors arrays, computing zeroth scattering 
      Intensity=zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
      for e in range(NEnergy):
            # print(e,log(Planck(x[e-1])/Planck(x[e]))/log(x[e]/x[e-1]),'           ',log(x[e])/log(10),'  ')
            for d in range(NZenith):
                  Intensity[e,d,0]=Planck(x[e])*(1 + 2.06*mu[d])
                  Intensity[e,d,1]=Intensity[e,d,0]*0.1171*(mu[d] - 1.)/(1. + 3.582*abs(mu[d]))
      s,gn=pr(s,gn,'I0')
      # exit()



# In[ ]:




if Spectrum=='FromFile' :
      inI = open(AtmName+'I.bin')
      inx = open(AtmName+'x.bin')
      inm = open(AtmName+'m.bin')
      x=fromfile(inx)
      mu=fromfile(inm)
      NEnergy=len(x)
      NZenith=len(mu)
      Intensity=fromfile(inI).reshape((NEnergy,NZenith,2))
      s,gn=pr(s,gn,'I is read')
else: 
      outI = open(AtmName+'I.bin','w')
      outx = open(AtmName+'x.bin','w')
      outm = open(AtmName+'m.bin','w')
      Intensity.tofile(outI,format="%e")
      x.tofile(outx,format="%e")
      mu.tofile(outm,format="%e")



# In[28]:



if plotAtm: # plot Everything and save All pics and tabs
      Iin=Planck # Delta # initial photon distribution 
      Stokes_defined=Spectrum=='Compton' or Spectrum=='Thomson'
      Sangles= range(NMu) if Stokes_defined else [] # list of the angle indeces to be plot a detailed figure
      Senergy=[int(0.5+n*log(1+4*Theta)/x_weight) for n in range(NColors) ] if Stokes_defined else [] 
      print(Senergy)
      outF = open(AtmName+'Fx.dat','w')
      outp = open(AtmName+'Pd.dat','w')
      frmt=lambda val, list: '{:>8}'.format(val)+': '+' '.join('{:15.5e}'.format(v) for v in list) +'\n'
      outp.write(frmt('Energies',x) )      
      outF.write(frmt('Energies',x) )       
      
      labelsize=20
      fontsize=25
      figA=figure(1,figsize=(10,11))
      figA.clear()
      figA.suptitle(r'$\tau_T=$'+str(tau_T)+
                    r'$,\,T={:5.1f}keV$'.format(T*evere/1e3)+
                    r'$,\,\Theta=$'+str(Theta),fontsize=fontsize)  
      plotAF=figA.add_subplot(2,1,1,xscale='log',yscale='log') 
      plotAp=figA.add_subplot(2,1,2,xscale='log')      
      
      xIinx=[(Iin(x[e])*x[e]) for e in range(NEnergy)]
      # plotAF.set_xlim([x[0],x[-1]])
      plotAF.set_xlim([1e-3,1.])

      plotAF.set_ylim([1e-3,1e-2])
      plotAF.set_ylabel(r'$xI_x(\tau_T,x)$',fontsize=fontsize)
      plotAF.tick_params(axis='both', which='major', labelsize=labelsize)
      plotAF.plot(x,xIinx,'k-.')
      
      plotAp.set_xlim([x[0],x[-1]])
      plotAp.tick_params(axis='both', which='major', labelsize=labelsize)
      plotAp.set_xlabel(r'$x\,[m_ec^2]$',fontsize=fontsize)
      plotAp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
      plotAp.plot(x,[.0]*NEnergy,'-.',color='xkcd:brown')  

      for d in range(NMu):
            d1=d+NMu
            z=str(int(arccos(mu[d1])*180/pi))
            xFx=(Intensity[:,d1,0]*x[:]) 
            plotAF.plot(x,xFx,colors[(d*NColors)//NMu])
            outF.write( frmt(z+'deg',xFx) )
            p=(Intensity[:,d1,1]/Intensity[:,d1,0]*1e2) 
            plotAp.plot(x,p,colors[(d*NColors)//NMu])
            outp.write( frmt(z+'deg',p) )
            if d in Sangles: # Specific angle 
                  figS=figure(2,figsize=(10,13))
                  figS.suptitle(r'$\tau_T=$'+str(tau_T)+
                                r'$,\,T={:5.1f}keV$'.format(T*evere*1e-3)+
                                r'$,\,\Theta=$'+str(Theta)+
                                r'$,\,\mu={:5.3f}$'.format(mu[d1])+
                                r'$\,(z\approx$'+z+
                                r'$^{\circ})$',fontsize=fontsize)  
                  plotSF=figS.add_subplot(3,1,1,xscale='log',yscale='log') 
                  plotSc=figS.add_subplot(3,1,2,xscale='log') 
                  plotSp=figS.add_subplot(3,1,3,xscale='log')      
      
                  plotSF.set_ylabel(r'$xI_x(\tau_T,x)$',fontsize=fontsize)
                  plotSF.tick_params(axis='both', which='major', labelsize=labelsize)
                  plotSF.set_xlim([x[0],x[-1]])
                  plotSF.set_ylim([1e-6,1e-1])
                  plotSF.plot(x,xFx,color='k',linewidth=2.1)
                  plotSF.plot(x,xIinx,'-.',color='xkcd:brown')
                  outF.write(frmt('Sc.N.0',xFx) )      
                  
                  plotSc.set_xlim([x[0],x[-1]])
                  plotSc.tick_params(axis='both', which='major', labelsize=labelsize)
                  plotSc.set_ylabel(r'$c\,[ \% ]$',fontsize=fontsize)
                  
                  plotSp.set_xlim([x[0],x[-1]])
                  plotSp.tick_params(axis='both', which='major', labelsize=labelsize)
                  plotSp.set_xlabel(r'$x\,[m_ec^2]$',fontsize=fontsize)
                  plotSp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
                  plotSp.plot(x,p,color='k',linewidth=2.1)
                  
                  for k in range(ScatterNum+1):
                        xFx=(Stokes_out[k,:,d1,0]*x[:]) 
                        contribution=(Stokes_out[k,:,d1,0]/Intensity[:,d1,0]*1e2) 
                        p=[.0]*NEnergy  if k==0 else (Stokes_out[k,:,d1,1]/Stokes_out[k,:,d1,0]*1e2) 
                        outF.write( frmt('Sc.N.'+str(k),xFx) )
                        outp.write( frmt('Sc.N.'+str(k),p) )
                        plotSF.plot(x,xFx,'--',color=colors[(k*NColors)//(ScatterNum+1)])
                        plotSc.plot(x,contribution,'--',color=colors[(k*NColors)//(ScatterNum+1)])      
                        plotSp.plot(x,p,'--',color=colors[(k*NColors)//(ScatterNum+1)])
                  figS.savefig(AtmName+'z'+z+'.pdf')
                  figS.clf()
      
      figA.savefig(AtmName+'zAll.pdf')
      figA.clf()
      figA.suptitle(r'$\tau_T=$'+str(tau_T)+
                    r'$,\,T={:5.1f}keV$'.format(T*evere/1e3)+
                    r'$,\,\Theta=$'+str(Theta),fontsize=fontsize)  
      plotAF=figA.add_subplot(2,1,1,xscale='linear',yscale='linear') 
      plotAp=figA.add_subplot(2,1,2,xscale='linear')      
      
      plotAF.set_xlim(0,1)
      # plotAF.set_ylim()
      plotAF.set_ylabel(r'$I_x(\tau_T,\mu)/I_x(\tau_T,1)$',fontsize=fontsize)
      plotAF.tick_params(axis='both', which='major', labelsize=labelsize)
      
      plotAp.set_xlim(0,1)
      plotAp.tick_params(axis='both', which='major', labelsize=labelsize)
      plotAp.set_xlabel(r'$\mu$',fontsize=fontsize)
      plotAp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
      plotAp.plot(mu[NMu:NZenith],[.0]*NMu,'-.',color='xkcd:brown')  



      # show() 
      outp.write('\n\n' )      
      outF.write('\n\n' )
      outp.write(frmt('Cosines',mu[NMu:]) )      
      outF.write(frmt('Cosines',mu[NMu:]) ) 
      outp.write(frmt('Angles',arccos(mu[NMu:])*180/pi) )
      outF.write(frmt('Angles',arccos(mu[NMu:])*180/pi) )
      k=-1     
      for e in range(NEnergy):
            Esp="{:<8.2e}".format(x[e])
            Inorm=(Intensity[e,NMu:,0]/Intensity[e,-1,0]) 
            p=(Intensity[e,NMu:,1]/Intensity[e,NMu:,0]*1e2) 
            outF.write( frmt(Esp,Inorm) )
            outp.write( frmt(Esp,p) )
            if e in Senergy:
                  k=k+1
                  # plotAp.plot(mu[NMu:NZenith],p,'-.',color='xkcd:black')  
                  # plotAF.plot(mu[NMu:NZenith],Inorm,'-.',color='xkcd:black')  
                  Inorm=(Stokes_out[k,e,NMu:,0]/Stokes_out[k,e,-1,0]) 
                  p=[.0]*NMu  if k==0 else (Stokes_out[k,e,NMu:,1]/Stokes_out[k,e,NMu:,0]*1e2) 
                  plotAp.plot(mu[NMu:],p,'-',color=colors[k])  
                  plotAF.plot(mu[NMu:],Inorm,'-',color=colors[k])
      figA.savefig(AtmName+'In.pdf')

      outF.close()
      outp.close()

      s,gn=pr(s,gn,'plap')



# In[30]:



if computePulse:


      NPhase = 120 # Number of equidistant phase points
      NBend= 20 # Number of knots in light bending integrations
      NAlpha= 1000 # 10000 # Number of psi/aplha grid points 
      IntBend = leggauss(NBend)

      phi,phi_weight=linspace(0,2*pi,num=NPhase,endpoint=False,retstep=True)
      nu=600 # star rotation frequency in Hz
      M=1.4 # star mass in solar masses
      R_g=M*2.95 # gravitational Schwarzschild radius
      R_e=12.0 # equatorial radius of the star in kilometers
      
      if oblateness=='AlGendy': # from AlGendy et. al. (2014)
            Omega_bar=2*pi*nu*sqrt(2*R_e**3/R_g)/c
            print('_O_^2',Omega_bar**2,'_O_',Omega_bar)
            flattening=(0.788-0.515*R_g/R_e)*Omega_bar**2 
            print(R_e*(1-flattening))
      elif oblateness=='Sphere':
            flattening=0.0
      else:
            flattening=oblateness
      # exit()
      print(flattening)

      outParams = open(PulsName+'.dat','w')
      outParams.write(AtmName+'.dat is the name of file with some corresponding atmosphere model\n')
      outParams.write(str(R_e)+' = equatorial radius R_e\n')
      outParams.write(str(M)+' = star mass M in solar masses\n')
      outParams.write(str(nu)+' = star rotation frequency nu in Hz\n')
      outParams.write(str(NPhase)+' '+str(NBend)+' '+str(NAlpha)+' = NPhase NBend NAlpha parameters\n')

      

      def Beloborodov(cos_psi):
            """Beloborodov's approximation for cos_alpha(cos_psi) light bending function
            takes the cos psi 
            returns the cos alpha and its derivative
            """
            return 1. + (cos_psi - 1.)/redshift**2 ,1./redshift**2

      def Schwarzschild(R,alpha):
            """Schwarzschild exact relation between the \psi and \\alpha angles, where
            \\alpha is the angle between radius vector of the spot and the direction of the outgoing photon near the surface
            and \psi is the angle between normal and light propagation at the limit of infinite distance.
            For given distance from the mass center and the emission angle \\alpha 
            this function returns two numbers: 
                  the corresponding angle \psi 
                  and the time lag over against the fotons emited with zero impact parameter at the radius.
            """
            kx,wx=IntBend
            eps=(1+kx[0])/4e2
            u=R_g/R 
            b=sin(alpha)/sqrt(1-u)*R # impact parameter
            if 2*alpha>pi+eps:
                  cos_3eta=sqrt(27)*R_g/2/b
                  if cos_3eta > 1:
                        return pi+2*eps,0 # the timelag 
                  closest_approach=-2*b/sqrt(3)*cos(arccos(cos_3eta)/3 + 2*pi/3)
                  psi_max, lag_max= Schwarzschild(closest_approach,pi/2.)
                  psi_min, lag_min= Schwarzschild(R,pi-alpha)
                  psi=2*psi_max - psi_min    
                  lag=2*lag_max - lag_min # + 2*(R - closest_approach + R_g*log((R - R_g)/(closest_approach - R_g)))/c 
                  if psi>pi:
                        return pi+eps,lag
            else:
                  psi=0
                  lag=(R_e - R + R_g*log( (R_e - R_g)/(R - R_g) ) )/c
                  for i in range(NBend):
                        ex=(kx[i]+1)/2
                        q=(2. - ex*ex - u*(1 - ex*ex)**2/(1 - u))*sin(alpha)**2
                        sr=sqrt(cos(alpha)**2+ex*ex*q)
                        if  2*alpha>pi-eps:
                              dpsi=b/R/sqrt(q)*wx[i] #*2/2
                        else:
                              dpsi=ex*b/R/sr*wx[i] #*2/2
                        dlag=dpsi*b/c/(1+sr)*wx[i] #*2/2
                        psi+= dpsi
                        lag+= dlag
            return psi,lag
      # flattening=0

      NRadius=2 + int(flattening*R_e/1e-1)
      print(NRadius)
      r, dr = linspace(R_e*(1 - flattening),R_e,num=NRadius,retstep=True)
      alpha, dalpha = linspace(0,arccos(-1/sqrt(2*r[0]/R_g/3)),NAlpha,retstep=True)
      psi=zeros((NRadius,NAlpha))
      dt=zeros((NRadius,NAlpha))
      for d in range(NRadius):
            print(d)
            for a in range(NAlpha):
                  psi[d,a],dt[d,a]=Schwarzschild(r[d],alpha[a])

      s,gn=pr(s,gn,'psidt')      
   

      Flux=zeros((NPhase,NEnergy,3))

      i=pi*4/18    # line of sight colatitude

      NSpots= 2*4 # * somewhat
      # theta = [pi/3,2*pi/3] # spot colatitude
      # l=[0,pi] # spot longitude
      # dS=[1,1] # some arbitrary units
      little=1e-2
      dS=[1]*NSpots
      l=[0,0,little,little,pi,pi,pi+little,pi+little]
      theta=[pi/3,pi/3+little,pi/3,pi/3+little,2*pi/3,2*pi/3+little,2*pi/3,2*pi/3+little]

      outParams.write(str(i)+' = sight colatitude i\n')
      outParams.write(str(theta)+' = spot colatitudes theta\n')
      outParams.write(str(l)+' = spot longitudes l\n')
      
      sin_i=sin(i)
      cos_i=cos(i)

      BoloFlux=zeros((NPhase,3))
      logIntensity=zeros((NEnergy,NZenith,4))
      logIntensity[:,:,:2] = log(absolute(Intensity))
      logIntensity[:,:,2:] = sign(Intensity)


      s,gn=pr(s,gn,'second precomp')
      for p in range(NSpots):
            sin_theta=sin(theta[p])
            cos_theta=cos(theta[p])
            print(flattening)
            R=R_e*(1 - flattening*cos_theta**2) 
            dR=R_e*flattening*cos_theta*sin_theta # dR / d\theta


            r1=bisect(r[1:-1],R) 
            r2=r1 + 1
            dr1=(R - r[r1])/dr
            dr2=(r[r2] - R)/dr
            # dr1,dr2=.5,.5


            redshift=1.0/sqrt(1.0 - R_g/R) # 1/sqrt(1-R_g/R) = 1+ z = redshift
            f=redshift/R*dR
            sin_gamma=f/sqrt(1 + f**2) # angle gamma is positive towards the north pole 
            cos_gamma=1.0/sqrt(1 + f**2)
            beta=2*pi*nu*R*redshift*sin_theta/c
            # beta=0
            # print(beta,nu,R,redshift,sin_theta)
            Gamma=1.0/sqrt(1.0 - beta**2)
            # exit()

            # alpha=linspace(0,arccos(-1/sqrt(2*R/R_g/3)),NAlpha)
            # psi=zeros(NAlpha)
            # dt=zeros(NAlpha)
            # for a in range(NAlpha):
            #       psi[a],dt[a]=Schwarzschild(R,alpha[a])
                  
            for t in range(NPhase):
                  if True: # find mu
                        sin_phi=sin(phi[t]+l[p])
                        cos_phi=cos(phi[t]+l[p])
                        cos_psi=cos_i*cos_theta + sin_i*sin_theta*cos_phi
                        sin_psi=sqrt(1. - cos_psi**2)
                        

                        psi0=arccos(cos_psi) 
                        a1=bisect(psi[r1], psi0)
                        a2=bisect(psi[r2], psi0)
                        psi1=psi[r1, a1]
                        psi2=psi[r2, a2]
                        dpsi1=psi1 - psi[r1, a1 - 1]
                        dpsi2=psi2 - psi[r2, a2 - 1]
                        dpsi=dpsi1*dr2 + dpsi2*dr1
                        dalpha1 = dalpha*(psi1 - psi0)/dpsi1
                        dalpha2 = dalpha*(psi2 - psi0)/dpsi2
                        alpha1=alpha[a1] - dalpha1
                        alpha2=alpha[a2] - dalpha2
                        
                        cos_alpha = cos(alpha2*dr1 + alpha1*dr2) # linear interpolation of alpha(psi)
                        sin_alpha = sqrt(1. - cos_alpha**2)
                        sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-4 else 1./redshift
                        dcos_alpha=sin_alpha_over_sin_psi *dalpha/dpsi # d cos\alpha \over d \cos \psi
                       
                        # cos_alpha, dcos_alpha=Beloborodov(cos_psi) # insert exact formula here
                        # sin_alpha = sqrt(1. - cos_alpha**2)
                        # sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-4 else 1./redshift
                        dt1=dt[r1,a1 - 1]*dalpha1/dalpha + dt[r1,a1]*(1. - dalpha1/dalpha)
                        dt2=dt[r2,a2 - 1]*dalpha2/dalpha + dt[r2,a2]*(1. - dalpha2/dalpha)
                        
                        dphi=(dt1*dr2 + dt2*dr1)*2*pi*nu # \delta\phi = \phi_{obs} - \phi 
                        # print(psi0,cos_alpha,dcos_alpha,dphi)
                        
                        # dphi=0

                        i0=dphi/phi_weight 
                        di1=i0 - floor(i0)
                        di2=floor(i0) + 1 - i0
                        t2=int((t + i0)%NPhase)
                        t1=t2 - 1

                        cos_xi = - sin_alpha_over_sin_psi*sin_i*sin_phi
                        delta = 1./Gamma/(1.-beta*cos_xi)
                        cos_sigma = cos_gamma*cos_alpha + sin_alpha_over_sin_psi*sin_gamma*(cos_i*sin_theta - sin_i*cos_theta*cos_phi)

                        sin_sigma = sqrt(1. - cos_sigma)
                        mu0=delta*cos_sigma # cos(sigma')
                        Omega=dS[p]*mu0*redshift**2*dcos_alpha #/9*16   
                        # Omegaarray[t]=max(Omega,0)
                        # print(t,' : \t',mu0,' \t ',dcos_alpha,'\t',dphi,cos_alpha,cos_psi,Omega)
                        if mu0<0: # this only for speeding up. the backwards intensity is usually zero
                              continue 


                  if True: # find chi
                        sin_chi_0= - sin_theta*sin_phi # times sin psi
                        cos_chi_0=sin_i*cos_theta - sin_theta*cos_i*cos_phi # times sin psi 
                        chi_0=arctan2(sin_chi_0,cos_chi_0)
                        
                        sin_chi_1=sin_gamma*sin_i*sin_phi*sin_alpha_over_sin_psi #times sin alpha sin sigma 
                        cos_chi_1=cos_gamma - cos_alpha*cos_sigma  #times sin alpha sin sigma 
                        chi_1=arctan2(sin_chi_1,cos_chi_1)
                        
                        sin_lambda=sin_theta*cos_gamma - sin_gamma*cos_theta
                        cos_lambda=cos_theta*cos_gamma + sin_theta*sin_gamma
                        cos_eps = sin_alpha_over_sin_psi*(cos_i*sin_lambda - sin_i*cos_lambda*cos_phi + cos_psi*sin_gamma) - cos_alpha*sin_gamma
                        # alt_cos_eps=(cos_sigma*cos_gamma - cos_alpha)/sin_gamma # legit! thanks God I checked it!
                        sin_chi_prime=cos_eps*mu0*Gamma*beta # times something
                        cos_chi_prime=1. - cos_sigma**2 /(1. - beta*cos_xi) # times the samething
                        chi_prime=arctan2(sin_chi_prime,cos_chi_prime)   

                        chi=chi_0 + chi_1 + chi_prime

                        sin_2chi=sin(2*chi)
                        cos_2chi=cos(2*chi)
                        # print(chi_prime,' \t',cos_chi_prime )


                  d2=bisect(mu[:-1],mu0)
                  d1=d2-1
                  # print(mu0,mu[d2],mu[d1],d2,d1,'       ')
                  mu1,mu2=mu[d1],mu[d2]
                  dmu, dmu1, dmu2 = mu2 - mu1, mu0 - mu1, mu2 - mu0
                  shift=delta/redshift


                  for e in range(NEnergy): 
                        x0=x[e]/shift
                        e1=bisect(x[1:-1],x0) # not the fastest way? anybody cares? ## seems, that light bending is more time consuming anyways
                        e2=e1+1
                        
                        x1, x2 = x[e1], x[e2]
                        dx, dx1, dx2 = x2 - x1, x0 - x1, x2 - x0
                        logIQ = (
                              dx2*dmu2*logIntensity[e1, d1] + 
                              dx2*dmu1*logIntensity[e1, d2] +
                              dx1*dmu2*logIntensity[e2, d1] +
                              dx1*dmu1*logIntensity[e2, d2] # bilinear interpolation of the Stokes parameters
                        )/dx/dmu 
                        I,Q=exp(logIQ[:2])* shift**3 * Omega* logIQ[2:]
                        # if any(absolute(logIQ[2:])<.9) :
                        #       I,Q = (
                        #       dx2*dmu2*Intensity[e1, d1] + 
                        #       dx2*dmu1*Intensity[e1, d2] +
                        #       dx1*dmu2*Intensity[e2, d1] +
                        #       dx1*dmu1*Intensity[e2, d2] # bilinear interpolation of the Stokes parameters
                        #       )/dx/dmu * shift**3 * Omega
                              
                        # print(dx,x_weight)
                        if I<0: ############
                              print('never')
                              # e_max=min(e_max,e)
                        F=array([I, Q*cos_2chi, Q*sin_2chi])
                        Flux[t2,e] += F*di1 # phase linear interpotalion
                        Flux[t1,e] += F*di2 

      s,gn=pr(s,gn,'curves done ')
      # print('trtrt',e_max)


# In[31]:



if plotPulse :
      outF = open(PulsName + 'F.bin','w')
      outf = open(PulsName + 'f.bin','w')
      Flux.tofile(outF,format="%e")
      phi.tofile(outf,format="%e")

      labelsize=20
      fontsize=25
      
      for e in range(0,NEnergy,2): # too many pictures
            phase=list(phi/2/pi)+[1.]
            I=zeros(NPhase+1)
            Q=zeros(NPhase+1)
            U=zeros(NPhase+1)
            for t in range(NPhase+1):
                  I[t],Q[t],U[t]=Flux[t-1,e]*x[e]

            p=sqrt(Q**2+U**2)/I*100
            PA=arctan2(-U,-Q)*90/pi+90
            
            figA=figure(2,figsize=(10,14))
            figA.clear()
            figA.suptitle(r'$\nu={:5.0f}Hz$'.format(nu)+
                          r'$,\,R_e={:5.1f}km$'.format(R_e)+
                          r'$,\,M=$'+str(M)+r'$M_{\odot}$'+',\n'+
                          r'$\,\theta={:5.1f}\degree$'.format(theta[0]*180/pi)+
                          r'$,\,i={:5.1f}\degree$'.format(i*180/pi)+',\n'+
                          r'$\,\lg(x)={:6.2f}$'.format(log(x[e])/log(1e1)),fontsize=fontsize)  
            plotAF=figA.add_subplot(3,1,1,yscale='linear') 
            plotAp=figA.add_subplot(3,1,2)      #
            plotAc=figA.add_subplot(3,1,3)      #

            # plotAF.set_xlim([x[0],x[-1]])
            plotAF.set_xlim(0,1)
            
            # plotAF.locator_params(axis='y', nbins=10)
            plotAF.set_ylabel(r'$F_x(\varphi)/F_x^{max}$',fontsize=fontsize)
            plotAF.tick_params(axis='both', which='major', labelsize=labelsize)
            # plotAF.plot(x,xIinx,'k-.')
            
            # plotAF.set_xlim([x[0],x[-1]])
            plotAp.set_xlim(0,1)
            plotAp.tick_params(axis='both', which='major', labelsize=labelsize)
            plotAp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
            plotAc.set_xlim(0,1)
            plotAc.set_ylim(0,180)
            plotAc.set_yticks([0,30,60,90,120,150,180])
            plotAc.tick_params(axis='both', which='major', labelsize=labelsize)
            plotAc.set_ylabel(r'$\chi\,[\degree]$',fontsize=fontsize)
            plotAc.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)
            
            col=colors[(e*NColors)//NEnergy]
            plotAF.plot(phase,I/I.max(),color=col)
            plotAp.plot(phase,p,color=col)
            plotAc.plot(phase,PA,color=col)

            figA.savefig(PulsName+'Ff{:03d}.pdf'.format(e))
            figA.clf()
            # figA.close()
      # show()

s,gn=pr(s,gn,'phase pics drawn ')      


# In[ ]:


print('end')

