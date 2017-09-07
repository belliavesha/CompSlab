import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import scipy.special as special
from scipy.interpolate import interp1d,interp2d
from scipy.integrate import quad
from numpy import exp, log
import matplotlib.pyplot as plt
import itertools as it

colors='rygcbm' # Rainbow
interp1dKind='cubic' # interpKind for interp1d
interp2dKind='cubic' # interpKind for interp2d 

# h = 6.626e-34
# c = 3.0e+8
# k = 1.38e-23
# keVinK=8.61732814974056e-8
# keVinHz=4.135665538536e-18
# Temperature = 10
# def planck(logE, T):
#     """Log E in keV and temperature in keV.
#     returns B(nu,T), where nu=E/h
        
#     """
#     E=10**logE
#     a = 2.0*(E/c/h)**2
#     intensity = (a*E/ (exp(E/T) - 1.0) if E/T>1e-5 else a*T) 
#     #print E,logE 
#     return intensity*4.11274466039867e-48#in joules
nu=1
  
I0=lambda nu : 1 # planck(nu,Temperature)
tI=lambda tau,mu:  (exp(-tau/mu)*I0(nu)) if mu>0 else 0 # zeroth scattering 
                                                
uA=lambda t,m: (1-m*m)*(Il[-1](t,m)+Il[-1](t,-m))*3/4 
uB=lambda t,m: m*m*(Ir[-1](t,m)+Ir[-1](t,-m))*3/8
uC=lambda t,m: (Ir[-1](t,m)+Ir[-1](t,-m))*3/8 # integrand functions
                         
uIl=lambda t,m : ( (1-m*m)*A[-1](t) + m*m*( B[-1](t)+C[-1](t) ) )
uIr=lambda t,m : ( B[-1](t)+C[-1](t) )

tA=lambda t: quad( (lambda mu: uA(t,mu)) , 0 , 1)[0] 
tB=lambda t: quad( (lambda mu: uB(t,mu)) , 0 , 1)[0] 
tC=lambda t: quad( (lambda mu: uC(t,mu)) , 0 , 1)[0] #integrations

tIl=np.vectorize(lambda t,m : uIl(t,          m) if abs(m)<meps \
       else quad( lambda t1 : uIl(t+m*log(t1),m), exp(((0 if m>0 else tau0)-t)/m) ,1 ) [0])  
tIr=np.vectorize(lambda t,m : uIr(t,          m) if abs(m)<meps \
       else quad( lambda t1 : uIr(t+m*log(t1),m), exp(((0 if m>0 else tau0)-t)/m) ,1 ) [0])  

muN = 40# Number of mu points
mur = np.linspace(0,1,num=muN) # mu range for plotting (only outward direction)    
mu =np.linspace(-1,1,num=muN*2+1) #  cosines range of emergent ray angles 
meps=mu[muN+1]/3   # for (mu < meps => m==0 ) condition
  
numberOfScatterings = 5 # and number of lines on the figure
tau0range=[.5,1.0,3.]  # different tau0 parameters
 # number of iterations (scatterings)
for it in range(len(tau0range)):
  tau0 = tau0range[it] # optical depth of the atmosphere
  tauN =8 # Number of tau points
  tau=np.linspace(0,tau0,num=tauN)  # optical depths range
  tau2d, mu2d = np.meshgrid(tau,mu)
  Il=[tI] #intensities polarized in meridional plane 
  Ir=[tI] # and the one that perpendicular to it

  A=[] # initializations
  C=[]
  B=[]                     

  for scatteringNumber in range(numberOfScatterings):
    print it, scatteringNumber
    f = open("tau"+str(tau0)+"ns"+str(scatteringNumber)+".dat","w")
    f.write("#muN= "+str(muN)+" ; mu:       Intensity:   polarization degree:\n")
    
    A.append(interp1d(tau,map(tA,tau), kind=interp1dKind)) # redefinition ofa integrations by interolations
    B.append(interp1d(tau,map(tB,tau), kind=interp1dKind)) # do.
    C.append(interp1d(tau,map(tC,tau), kind=interp1dKind)) # do.
    
    Il.append(interp2d(tau2d,mu2d,tIl(tau2d,mu2d),kind=interp2dKind)) # ditto
    Ir.append(interp2d(tau2d,mu2d,tIr(tau2d,mu2d),kind=interp2dKind)) # do.
    
    N=Il[-1](tau0,1)+Ir[-1](tau0,1) 
    I=Il[-1](tau0,mur)+Ir[-1](tau0,mur)
    lI=Il[-1](tau0,mur)/N
    rI=Ir[-1](tau0,mur)/N
    Q=(Il[-1](tau0,mur)-Ir[-1](tau0,mur))/I

    for l in range(muN):
        f.write(\
            str(mu[l]).ljust(20) + \
            str(I[l][0]).ljust(20) + \
            str(Q[l][0]).ljust(20)+"\n")
   
    plt.figure(1)

    plt.subplot(2,3,it+1)
    plt.title(r'$\tau_0=$'+str(tau0))
    if it==0: 
        plt.ylabel(r'${I(\tau_0,\mu)}/{I(\tau_0,1)}$')
    plt.plot(mur, I/N, colors[scatteringNumber])

    plt.subplot(2,3,it+4)
    if it==0: 
        plt.ylabel(r'$p\,[ \% ]$')
    plt.xlabel(r'$\mu$')
    plt.plot(mur, 1e2*Q, colors[scatteringNumber])
    f.close()

caption="""The scattering numbers are indicated by colors:
1st --- red; 2nd --- yellow; 3rd --- green; 4th --- cyan; 5th --- blue."""  #6th --- purple
#plt.figtext(.1, .1,caption)   
print caption

plt.show()
    
print 'end'    
    


