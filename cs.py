# timer:
from time import time
time0=time()
def pr(s,gn,message = ''): 
      print gn,' :'+message+': ',time()-s
      return time(),gn+1
s,gn=pr(time0,0)

#import:
from numpy import linspace, logspace, empty, zeros,ones
from numpy import pi, exp, log, sqrt, sin, cos
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.special import kn
from matplotlib.pyplot import *
#import numpy as np
#import matplotlib.pyplot as plt
s,gn=pr(s,gn, 'importing')

colors='rygcbm' # Rainbow

#physical constants:
evere=.5109989e6 # electron volts in elecron rest energy 
incgs=3.43670379e30 # 2 m_e^4 c^6 / h^3

#precomputations :
ScatterNum = 10 # total number of scatterings
NGamma=10 # number of Lorenz factor points (\gamma)
NAzimuth=10 # (*) numbers of azimuth angles (\phi)
NDirection = 14 # (*) number of propagation angle cosines (\mu) 
NEnergy = 50 # number of energy points (x)
NDepth = 10 # number of optical depth levels (\tau)
      # (*) Notice that if some numer of points is odd then there may be a zero angle or a cosine which equals to 1
      # That may lead to Zero Division Error in these functions, be careful and cautious 
tau_T= .5 # Thomson optical depth of thermalization 
x_l, x_u = -5 , .1 # lower and upper bounds of the log_10 energy span

IntGamma = laggauss(NGamma) # sample points and weights for computing thermal matrix
IntAzimuth = leggauss(NAzimuth) # sample points and weights for computing azimuth-averaged matrix
IntDirection = leggauss(NDirection) #  sample points and weights for computing the source function of one given incoming photon energy
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for computing the full source function 
IntDepth = linspace(0,tau_T,num=NDepth,retstep=True)  # sample points and weights for computing the intensities

def setTe(T):
      " Use this function to set electron gas Maxwellian temperature"
      K = kn(2,1./T) # second modified Bessel function of reversed dimensionless temperature       
      return (T, K) # Theta is the dimensionless electron gas temperature (Theta = k * T_e / m_e c^2)

Theta, K2Y = setTe(.5e5/evere) # it's about 0.1
T = 1e1/evere # photon black body temperature
s,gn=pr(s,gn,'precomps')


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

def sigma_cs(x): # if Theta isn't small one will need to compute the mean cross-section  based on electron momentum distribution
      """ This function compute the Compton scattering cross-section in electron rest frame 
      x is the energy of the scattering photon in units of electron rest energy
      this function approaches the mean compton cross-section when electron gas temperature is small
      """ 
      if x<.1:
            a,n,s=3./8.,0.,0.
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

      # print  gamma-gammaStar, gammaStar,gamma,Q,u,(u-Q)/(Q+u)

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

if False : #check cross section  
      cro = open('cros.dat','r')
      figure(10)
      x=[]
      cs=[]
      sgm=[]
      for n in range(85):
            line = cro.readline().split()
            de=lambda X : float (X[:X.find('D')]+'e'+X[X.find('D')+1:])
            x.append(de(line[3]))
            cs.append(de(line[5]))
            sgm.append(sigma_cs(10**(x[-1])))
            # print x, cs
      plot(x,cs)
      plot(x,sgm)
      show()
      exit()

if True : # Check spectra
      spe=open('tes1.spe','r')
      figure(1)
      for n in range(5):
            pra=spe.readline()
            # print pra
            x=[]
            # xFx1=[]
            # xFx2=[]
            # xFx3=[]
            # xFx4=[]
            xFx5=[]
            for e in range(85):
                  line=spe.readline().split()
                  x.append(float (line[0]))
                  # xFx1.append(line[1])
                  # xFx2.append(line[2])
                  # xFx3.append(line[3])
                  # xFx4.append(line[4])
                  xFx5.append(float (line[5]) * x[-1])
            # plot(x,xFx1,colors[1])
            # plot(x,xFx2,colors[2])
            # plot(x,xFx3,colors[3])
            # plot(x,xFx4,colors[4])
            plot(x,xFx5,colors[5])
            s,gn=pr(s,gn,pra[-2])
      xscale('log')
      yscale('log')
      # show()
      # exit()

if False: #draw Planck
      x=IntEnergy[0]
      print x
      print map(lambda e: (e*Planck(e)),x)
      
      plot(x,map(lambda e: (e*Planck(e)),x),'g')
      xscale('log')
      yscale('log')
      show()

if False: #draw delta
      x=IntEnergy[0]
      print map(lambda e: (e*Delta(e)),x)
      print x
      plot(x,map(lambda e: (e*Delta(e)),x),'g')
      xscale('log')
      yscale('log')
      show()
      exit()

if False : # compare the numbers to ones obtained from old fortran code
      from compare import Factor
      Factor(Compton_redistribution_aa)
      
if False : # checking symmetries
      from check import *
      print CheckAngularSymmetry(Compton_redistribution_aa,0.01,0.1,-0.4,0.5)
      s,gn=pr(s,gn,'ang-check')
      print CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.1,0.5,-0.5,Theta)
      s,gn=pr(s,gn,'freq-check')
      exit()

if True: # Computing Compton scattering cross section for all energies
      x,x_weight=IntEnergy
      sigma=map(sigma_cs,x)   
      # sigma = zeros(NEnergy) + 1

if True : # Computing redistribution matrices for all energies and angles 
      mu,mu_weight=IntDirection
      RedistributionMatrix = ones( (NEnergy,NEnergy,NDirection,NDirection,2,2) )
      for e in range(NEnergy): # x
            for e1 in range(e,NEnergy): # x1
                  for d in range(NDirection/2): # mu 
                        for d1 in range(d,NDirection/2): #mu1
                              md=NDirection-d-1 # -mu
                              md1=NDirection-d1-1 # -mu1
                              t=d1>d
                              f=e1>e

                              if 1: 
                                    r=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[d1])
                                    rm=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[md1])
                                    RedistributionMatrix[e][e1][d][d1]=r
                                    RedistributionMatrix[e][e1][md][md1]=r
                                    RedistributionMatrix[e][e1][d][md1]=rm
                                    RedistributionMatrix[e][e1][md][d1]=rm
                              if t: # angular symmetry
                                    rt=r # transponed
                                    rmt=rm # matrices
                                    rt[0][1],rt[1][0]=r[1][0],r[0][1]
                                    rmt[0][1],rmt[1][0]=rm[1][0],rm[0][1]
                                    RedistributionMatrix[e][e1][d1][d]=rt
                                    RedistributionMatrix[e][e1][md1][md]=rt
                                    RedistributionMatrix[e][e1][md1][d]=rmt
                                    RedistributionMatrix[e][e1][d1][md]=rmt
                              if f: # frequency symmetry
                                    m=exp((x[e]-x[e1])/Theta)
                                    rf=r*m  # when Maxwellian
                                    rmf=rm*m  # or Wein distributions
                                    RedistributionMatrix[e1][e][d][d1]=rf
                                    RedistributionMatrix[e1][e][md][md1]=rf
                                    RedistributionMatrix[e1][e][d][md1]=rmf
                                    RedistributionMatrix[e1][e][md][d1]=rmf
                              if t and f: # both symmeties 
                                    rtf=rt*m
                                    rmtf=rmt*m
                                    RedistributionMatrix[e1][e][d1][d]=rtf
                                    RedistributionMatrix[e1][e][md1][md]=rtf
                                    RedistributionMatrix[e1][e][md1][d]=rmtf
                                    RedistributionMatrix[e1][e][d1][md]=rmtf
      s,gn=pr(s,gn,'RMtable')

if True : # Initializing Stokes vectors arrays, computing zeroth scattering 
      Iin=Planck
      tau,tau_weight=IntDepth
      Source=zeros((ScatterNum,NDepth,NEnergy,NDirection,2)) # source function                 
      Stokes=zeros((ScatterNum,NDepth,NEnergy,NDirection,2)) # intensity Stokes vector
      Stokes_out=zeros((ScatterNum+1,NEnergy,NDirection,2)) # outgoing Stokes vector of each scattering
      Stokes_in=zeros((NDepth,NEnergy,NDirection,2)) # Stokes vector of the initial raiation (0th scattering) 
      Intensity=zeros((NEnergy,NDirection,2)) # total intensity of all scattering orders from the slab suface 
      for e in range(NEnergy):
            for d in range(NDirection):
                  for t in range(NDepth):
                        Stokes_in[t][e][d][0]=Iin(x[e])*exp(-tau[t]*sigma[e]/mu[d]) if mu[d]>0 else 0 
                        Stokes_in[t][e][d][1]=0
                  else:
                        Stokes_out[0][e][d][0]=Iin(x[e])*exp(-tau_T*sigma[e]/mu[d]) if mu[d]>0 else 0
                        Stokes_out[0][e][d][1]=0
      Intensity += Stokes_out[0]
      s,gn=pr(s,gn,'I0')

for k in range(ScatterNum): # do ScatterNum scattering iterations
      for t in range(NDepth): # S_k= R I_{k-1}
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
                                    if (r[0][0]-1.)*(r[0][0]-1.)<1e-10 : print '!!!', e,e1,d,d1  
                        Source[k][t][e][d]+=S      
      
      for t in range(NDepth):
            for e in range(NEnergy): # I_k= integral S_k
                  for d in range(NDirection): 
                        I=zeros(2)
                        for t1 in range(t+1) if mu[d]>0 else range (t,NDepth): #here is where zeros used to come from
                              S=Source[k][t1][e][d]
                              I+=tau_weight*S*exp(sigma[e]*(tau[t1]-tau[t])/mu[d])
                        Stokes[k][t][e][d]+=I/abs(mu[d]) #abs
                        if I[0]>100 : print 'wow!',t,e,d
      else:
            Stokes_out[k+1]+=Stokes[k][t]
            Intensity+=Stokes[k][t]    
      s,gn=pr(s,gn,'I'+str(1+k))

out = open('res','w')
d=NDirection-1#*3/5
out.write(str(mu[d])+str(list(x))+'\n')
#range(NDirection/2,NDirection):
print d, mu[d]

if True: # plotting intensity
      out.write('xFlux\n')
      figure(1)
      for k in range(ScatterNum+1):
            xFx=[(Stokes_out[k][e][d][0]*x[e]) for e in range(NEnergy)]
            plot(x,xFx,colors[(k*5)/ScatterNum])
            out.write(str(k)+' : '+str(mu[d])+str(xFx)+'\n')
      else: 
            xFx=[(Intensity[e][d][0]*x[e]) for e in range(NEnergy)]
            plot(x,xFx,'k')
            out.write('st : '+str(mu[d])+str(xFx)+'\n')
      yscale('log')
      xscale('log')
      s,gn=pr(s,gn,'plfl')

if False: # plotting polarization
      figure(2)
      out.write('ppc\n')
      p = lambda a :[1e2*a[e][d][1]/a[e][d][0] for e in range(NEnergy)]
      for k in range(1+ScatterNum):
            print k
            yr=p(Stokes_out[k])
            plot(x,yr,colors[(k*5)/ScatterNum])
            out.write(str(k)+' : '+str(mu[d])+str(yr)+'\n')
      else: 
            yr=p(Intensity) 
            plot(x,yr,'k')
            out.write('pt : '+str(mu[d])+str(yr)+'\n')
      xscale('log')
      s,gn=pr(s,gn,'plpl')

if False: # plotting intensity
      out.write('xFlux\n')
      figure(3)
      for d in range(NDirection/2,NDirection):
            xFx=[(Intensity[e][d][0]*x[e]) for e in range(NEnergy)]
            plot(x,xFx,colors [(5*(d-NDirection/2))*2/NDirection ])
            out.write(str(d)+' : '+str(mu[d])+str(xFx)+'\n')
      yscale('log')
      xscale('log')
      s,gn=pr(s,gn,'plfl')

if False: # plotting plarization
      figure(4)
      out.write('ppc\n')
      p = lambda a :[1e2*a[e][d][1]/a[e][d][0] for e in range(NEnergy)]
      for d in range(NDirection/2,NDirection):
            yr=p(Intensity) 
            plot(x,yr,colors [5*(d-NDirection/2)*2/NDirection ])
            out.write(str(d)+' : '+str(mu[d])+str(yr)+'\n')
      xscale('log')
      s,gn=pr(s,gn,'plpl')



# diff = (lambda a,b : (log(Intensity[a][d][0])-log(Intensity[b][d][0]))/(log(x[a])-log(x[b])))
# print diff(10,12),diff(12,14) ,diff(14,16),diff(16,18)  
show()  

print 'Total time : ', s-time0

