
# coding: utf-8

# In[1]:


# timer:
from time import time
time0=time()
def pr(s,gn,message = ''): 
      print(gn,' :'+message+': ',time()-s)
      return time(),gn+1
s,gn=pr(time0,0)


# In[2]:



#import:
from numpy import linspace, logspace, empty, zeros, ones
from numpy import pi, exp, log, sqrt, sin, cos, arccos
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.special import kn
from matplotlib.pyplot import *

#import numpy as np
#import matplotlib.pyplot as plt
s,gn=pr(s,gn, 'importing')


# In[3]:



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
# incgs=3.43670379e30 # 2 m_e^4 c^6 / h^3

# parameters: 
tau_T= 0.9# Thomson optical depth of thermalization 
x_l, x_u = -5.1 , 0 # lower and upper bounds of the log_10 energy span
Theta = 0.1 # dimensionless electron gas temperature (Theta = k T_e / m_e c^2) # it's about 0.1 
T = 10/evere #  dimensionless photon black body temperature T = k T_bb / m_e c^2
saveName='T10Thp1tp9' # the prefix for all result files related to the set of parameters
# maybe here better be some clever algorythm compiling T, Theta and tau_T in a proper and readable filename 

#precomputations :
ScatterNum = 10 #20 # total number of scatterings
NGamma= 10 # number of Lorenz factor points (\gamma)
NAzimuth= 10 # numbers of azimuth angles (\phi) [0,pi]
NEnergy = 21 # 91 # number of energy points (x)
NDepth = 9 # 41 # number of optical depth levels (\tau)
NMu = 6 # 15 # number of propagation zenith angle cosines (\mu) [0,1]
NZenith = 2*NMu # number of propagation zenith angles (z) [0,pi]

IntGamma = laggauss(NGamma) # sample points and weights for computing thermal matrix
IntAzimuth = leggauss(NAzimuth*2) # sample points and weights for computing azimuth-averaged matrix
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function
IntDepth = linspace(0,tau_T,num=NDepth,retstep=True)  # sample points and weights for integrations over the optical depth computing intencity 
IntZenith = leggauss(NZenith) #  sample points and weights for integrations over zenith angle in positive and negative directions together

K2Y = kn(2,1./Theta) # second modified Bessel function of reversed dimensionless temperature       

s,gn=pr(s,gn,'precomps')


# In[4]:



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
      and also the are non-zero elements of the matrix: R11,R12=R21,R22,R33 respectfully 
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
      
      az_c=cos(pi*phi)  # list of azimuth cosines
      az_s=sin(pi*phi)**2  # list of azimuth square sinuses
      sc_c=mu1*mu2-sqrt(eta1*eta2)*az_c # list of scattering angles' cosines
      sc_s=1. - sc_c**2 # list of scattering angles' squared sinuses
      cos2chi1 = 2.*(mu1*sc_c-mu2)*(mu1*sc_c-mu2)/eta1/sc_s-1.  # list[ cos( 2 \chi_1 ) ]
      cos2chi2 = 2.*(mu1-mu2*sc_c)*(mu1-mu2*sc_c)/eta2/sc_s-1.  # list[ cos( 2 \chi_2 ) ]
      sin2chiP = 4.*(mu1-mu2*sc_c)*(mu1*sc_c-mu2)*az_s/sc_s**2  # list[ sin( 2 \chi_1 )*sin( 2 \chi_2 ) ]

      R=zeros( (2,2,) )
      for i in range(NAzimuth*2):
            (C,I,Q,U)=Compton_redistribution(x1,x2,sc_c[i])
            R[0][0]+=C*pi*phi_weight[i]
            R[0][1]+=I*pi*cos2chi2[i]*phi_weight[i]
            R[1][0]+=I*pi*cos2chi1[i]*phi_weight[i]
            R[1][1]+=pi*(Q*cos2chi1[i]*cos2chi2[i]+U*sin2chiP[i])*phi_weight[i]         
      
      
      # print(x1,x2,mu1,mu2,R)
      return R*x1*x1/x2

s,gn=pr(s,gn,'funcs')


# In[5]:


if True : # checking symmetries #unnecessary 
      from check import *
      print('Angular symmetry: ',CheckAngularSymmetry(Compton_redistribution_aa,0.02,0.02,-0.4,0.5))
      s,gn=pr(s,gn,'ang-check')
      print('Energy symmetry: ',CheckFrequencySymmetry(Compton_redistribution_aa,0.01,0.02,0.5,0.5,Theta))
      s,gn=pr(s,gn,'freq-check')
      # exit()


# In[6]:





if True : # Computing redistribution matrices for all energies and angles 
      mu,mu_weight=IntZenith
      x,x_weight=IntEnergy
      sigma=zeros(NEnergy)
      RedistributionMatrix = ones( (NEnergy,NEnergy,NZenith,NZenith,2,2) )
      for e in range(NEnergy): # x [-\infty,\infty]
            for e1 in range(e,NEnergy): # x1 [x,\infty]
                  for d in range(NMu): # mu [-1,0]
                        for d1 in range(d,NMu): # mu1 [-1,mu]
                              md=NZenith-d-1 # -mu
                              md1=NZenith-d1-1 # -mu1
                              w=mu_weight[d1]*x_weight*mu_weight[d]
                              t=d1>d
                              f=e1>e

                              if 1: 
                                    r=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[d1])
                                    rm=Compton_redistribution_aa(x[e],x[e1],mu[d],mu[md1])
                                    sigma[e1]+=(r[0][0]+rm[0][0])*w
                                    RedistributionMatrix[e][e1][d][d1]=r
                                    RedistributionMatrix[e][e1][md][md1]=r
                                    RedistributionMatrix[e][e1][d][md1]=rm
                                    RedistributionMatrix[e][e1][md][d1]=rm
                              if t: # angular symmetry
                                    rt=r # transponed
                                    rmt=rm # matrices
                                    rt[0][1],rt[1][0]=r[1][0],r[0][1]
                                    rmt[0][1],rmt[1][0]=rm[1][0],rm[0][1]
                                    sigma[e1]+=(rt[0][0]+rmt[0][0])*w
                                    RedistributionMatrix[e][e1][d1][d]=rt
                                    RedistributionMatrix[e][e1][md1][md]=rt
                                    RedistributionMatrix[e][e1][md1][d]=rmt
                                    RedistributionMatrix[e][e1][d1][md]=rmt
                              if f: # frequency symmetry
                                    m=exp((x[e]-x[e1])/Theta)*x[e1]**3/x[e]**3
                                    rf=r*m  # when Maxwellian
                                    rmf=rm*m  # or Wein distributions
                                    sigma[e]+=(rf[0][0]+rmf[0][0])*w
                                    RedistributionMatrix[e1][e][d][d1]=rf
                                    RedistributionMatrix[e1][e][md][md1]=rf
                                    RedistributionMatrix[e1][e][d][md1]=rmf
                                    RedistributionMatrix[e1][e][md][d1]=rmf
                              if t and f: # both symmeties 
                                    rtf=rt*m
                                    rmtf=rmt*m
                                    sigma[e]+=(rtf[0][0]+rmtf[0][0])*w
                                    RedistributionMatrix[e1][e][d1][d]=rtf
                                    RedistributionMatrix[e1][e][md1][md]=rtf
                                    RedistributionMatrix[e1][e][md1][d]=rmtf
                                    RedistributionMatrix[e1][e][d1][md]=rmtf
      s,gn=pr(s,gn,'RMtable')



# In[11]:



if True : #check cross section  
#       cro = open('cros.dat','r')
      figure(1,figsize=(10,9))
      xscale('log')
#       fx=[0.0]
#       cs=[1.0]
#       sgm=[1.0]
#       de=lambda X : float (X[:X.find('D')]+'e'+X[X.find('D')+1:])
#       for n in range(85):
#             line = cro.readline().split()
#             fx.append(10**de(line[3]))
#             cs.append(de(line[5]))
#             sgm.append(sigma_cs(fx[-1]))
#             # print x, cs
#       cs.append(0)
#       fx.append(x[-1])
#       sgm.append(sigma_cs(x[-1]))
#       fcs=interp1d(fx,cs)
#       fcsx=fcs(x)
#       plot(fx,cs)
      plot(x,list(map(sigma_cs,x)),'b')
      plot(x,sigma,'r')
#       plot(x,fcsx,'k')
#       plot(fx,sgm)
#       plot(x,sigma2-fcsx)
#       print(sigma,list(map(sigma_cs,x)))
      savefig('compsigma.png')
      show()
      s,gn=pr(s,gn,'sigmaplot') 


# In[8]:



if True : # Initializing Stokes vectors arrays, computing zeroth scattering 
      Iin=Planck # Delta # initial photon distribution 
      tau,tau_weight=IntDepth
      Source=zeros((ScatterNum,NDepth,NEnergy,NZenith,2)) # source function                 
      Stokes=zeros((ScatterNum,NDepth,NEnergy,NZenith,2)) # intensity Stokes vector
      Stokes_out=zeros((ScatterNum+1,NEnergy,NZenith,2)) # outgoing Stokes vector of each scattering
      Stokes_in=zeros((NDepth,NEnergy,NZenith,2)) # Stokes vector of the initial raiation (0th scattering) 
      Intensity=zeros((NEnergy,NZenith,2)) # total intensity of all scattering orders from the slab suface 
      for e in range(NEnergy):
            for d in range(NZenith):
                  for t in range(NDepth):
                        Stokes_in[t][e][d][0]=Iin(x[e])*exp(-tau[t]*sigma[e]/mu[d]) if mu[d]>0 else 0 
                        Stokes_in[t][e][d][1]=0
                  else:
                        Stokes_out[0][e][d][0]=Iin(x[e])*exp(-tau_T*sigma[e]/mu[d]) if mu[d]>0 else 0
                        Stokes_out[0][e][d][1]=0
      Intensity += Stokes_out[0]
      s,gn=pr(s,gn,'I0')



# In[9]:



for k in range(ScatterNum): # do ScatterNum scattering iterations
      for t in range(NDepth): # S_k= R I_{k-1}
            for e in range(NEnergy):
                  for d in range(NZenith):
                        S=zeros(2)
                        for e1 in range(NEnergy):
                              for d1 in range(NZenith):
                                    w = mu_weight[d1]*x_weight # total weight
                                    r = RedistributionMatrix[e][e1][d][d1]  #  
                                    I = Stokes[k-1][t][e1][d1] if k>0 else Stokes_in[t][e1][d1]
                                    S[0]+= w*( I[0]*r[0][0] + I[1]*r[0][1] ) # 
                                    S[1]+= w*( I[0]*r[1][0] + I[1]*r[1][1] ) #
                        Source[k][t][e][d]+=S #     
      
      for t in range(NDepth):# I_k= integral S_k
            for e in range(NEnergy): 
                  for d in range(NZenith): 
                        I=zeros(2)
                        I+=tau_weight/2*Source[k][t][e][d]
                        for t1 in range(t) if mu[d]>0 else range (t+1,NDepth): #
                              S=Source[k][t1][e][d] #
                              I+=tau_weight*S*exp(sigma[e]*(tau[t1]-tau[t])/mu[d])
                        Stokes[k][t][e][d]+=I/abs(mu[d]) #abs
      else:
            Stokes_out[k+1]+=Stokes[k][t]
            Intensity+=Stokes[k][t]   
       
      s,gn=pr(s,gn,'I'+str(1+k))


# In[10]:


# ScatterNum=20
if True: # plot Everything and save All pics and tabs
      outF = open(saveName+'Fx.dat','w')
      outp = open(saveName+'Pd.dat','w')
      frmt=lambda val, list: '{:>8}'.format(val)+': '+' '.join('{:15.5e}'.format(v) for v in list) +'\n'
      outp.write(frmt('Energies',x) )      
      outF.write(frmt('Energies',x) )       
      
      labelsize=20
      fontsize=25
      figA=figure(1,figsize=(16,18))
      figA.suptitle(r'$\tau_T=$'+str(tau_T)+
                    r'$,\,T={:5.1f}eV$'.format(T*evere)+
                    r'$,\,\Theta=$'+str(Theta),fontsize=fontsize)  
      plotAF=figA.add_subplot(2,1,1,xscale='log',yscale='log') 
      plotAp=figA.add_subplot(2,1,2,xscale='log')      
      
      xIinx=[(Iin(x[e])*x[e]) for e in range(NEnergy)]
      plotAF.set_xlim([x[0],x[-1]])
      plotAF.set_ylim([1e-7,2e-4])
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
            xFx=[(Intensity[e][d1][0]*x[e]) for e in range(NEnergy)]
            plotAF.plot(x,xFx,colors[(d*NColors)//NMu])
            outF.write( frmt(z+'deg',xFx) )
            p=[(Intensity[e][d1][1]/Intensity[e][d1][0]*1e2) for e in range(NEnergy)]
            plotAp.plot(x,p,colors[(d*NColors)//NMu])
            outp.write( frmt(z+'deg',p) )
            if True: # Specific angle 
                  figS=figure(2+d,figsize=(16,21))
                  figS.suptitle(r'$\tau_T=$'+str(tau_T)+
                                r'$,\,T={:5.1f}eV$'.format(T*evere)+
                                r'$,\,\Theta=$'+str(Theta)+
                                r'$,\,\mu={:5.3f}$'.format(mu[d])+
                                r'$\,(z\approx$'+z+
                                r'$^{\circ})$',fontsize=fontsize)  
                  plotSF=figS.add_subplot(3,1,1,xscale='log',yscale='log') 
                  plotSc=figS.add_subplot(3,1,2,xscale='log') 
                  plotSp=figS.add_subplot(3,1,3,xscale='log')      
      
                  plotSF.set_ylabel(r'$xI_x(\tau_T,x)$',fontsize=fontsize)
                  plotSF.tick_params(axis='both', which='major', labelsize=labelsize)
                  plotSF.set_xlim([x[0],x[-1]])
                  plotSF.set_ylim([1e-17,3e-3])
                  plotSF.plot(x,xFx,'k')
                  plotSF.plot(x,xIinx,'-.',color='xkcd:brown')
                  outF.write(frmt('Sc.N.0',xFx) )      
                  
                  plotSc.set_xlim([x[0],x[-1]])
                  plotSc.tick_params(axis='both', which='major', labelsize=labelsize)
                  plotSc.set_ylabel(r'$c\,[ \% ]$',fontsize=fontsize)
                  
                  plotSp.set_xlim([x[0],x[-1]])
                  plotSp.tick_params(axis='both', which='major', labelsize=labelsize)
                  plotSp.set_xlabel(r'$x\,[m_ec^2]$',fontsize=fontsize)
                  plotSp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
                  plotSp.plot(x,p,'k')
                  
                  for k in range(ScatterNum+1):
                        xFx=[(Stokes_out[k][e][d1][0]*x[e]) for e in range(NEnergy)]
                        c=[(Stokes_out[k][e][d1][0]/Intensity[e][d1][0]*1e2) for e in range(NEnergy)]
                        p=[.0]*NEnergy  if k==0 else (
                            [(Stokes_out[k][e][d1][1]/Stokes_out[k][e][d1][0]*1e2) for e in range(NEnergy)])
                        outF.write( frmt('Sc.N.'+str(k),xFx) )
                        outp.write( frmt('Sc.N.'+str(k),p) )
                        plotSF.plot(x,xFx,'--',color=colors[(k*NColors)//(ScatterNum+1)])
                        plotSc.plot(x,c,'--',color=colors[(k*NColors)//(ScatterNum+1)])      
                        plotSp.plot(x,p,'--',color=colors[(k*NColors)//(ScatterNum+1)])
                  figS.savefig(saveName+'z'+z+'.eps')
                  figS.savefig(saveName+'z'+z+'.pdf')
      
      figA.savefig(saveName+'zAll.pdf')
      figA.savefig(saveName+'zAll.eps')  
      show() 
      outp.write('\n\n' )      
      outF.write('\n\n' )
      outp.write(frmt('Cosines',mu[NMu:]) )      
      outF.write(frmt('Cosines',mu[NMu:]) ) 
      outp.write(frmt('Angles',arccos(mu[NMu:])*180/pi) )
      outF.write(frmt('Angles',arccos(mu[NMu:])*180/pi) )      
      for e in range(NEnergy):
            Esp="{:<8.2e}".format(x[e])
            xFx=[(Intensity[e][d1][0]*x[e]) for d1 in range(NMu,NZenith)]
            p=[(Intensity[e][d1][1]/Intensity[e][d1][0]*1e2) for d1 in range(NMu,NZenith)]
            outF.write( frmt(Esp,xFx) )
            outp.write( frmt(Esp,p) )
      outF.close()
      outp.close()

      s,gn=pr(s,gn,'plap')
