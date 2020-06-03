#import:
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *

import numpy as np 
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil
# from numpy.polynomial.laguerre import laggauss
# from numpy.polynomial.legendre import leggauss
# from scipy.interpolate import interp1d
# from scipy.special import kn
# from bisect import bisect

colors=['xkcd:brownish red',
        'xkcd:red',
        # 'xkcd:orange',
        # 'xkcd:dark yellow',
        # 'xkcd:dark yellow green',
        'xkcd:deep green',
        # 'xkcd:dark cyan',
        # 'xkcd:blue',
        'xkcd:purple'   
]

#physical constants:
evere=.5109989e6 # electron volts in elecron rest energy 
G=13275412528e1 # G*M_sol in km^3/s^2 
c=299792458e-3 # speed of light in km/s

NPhi = 150#500#150#128*2#120#128#120 # Number of equidistant phase points
phi, phi_weight=linspace(0,360,num=NPhi,endpoint=False,retstep=True) #Size of spacing between samples = phi_weight
chi = zeros(NPhi)
chi_sph = zeros(NPhi)

nu=600#1#100#600#700#600 # star rotation frequency in Hz
M=1.4#1.6#1.4 # star mass in solar masses
R_g=M*2.95 # gravitational Schwarzschild radius
R_e=12.0 # equatorial radius of the star in kilometers      
incl = 5*pi/18 #pi*4/18#1.5    # line of sight colatitude #in radians
theta = [pi/2,3*pi/18,6*pi/18,pi/18]
NSpots = len(theta)



#rc("text", usetex=True)
figA = figure(figsize=(14,18), dpi=300) #8,6
#rc("font", family="serif")
#rc("font",serif="Times")
matplotlib.pyplot.figure(1)
lbfontsz = 25 
lwidth= 2.5
rc("xtick", labelsize=lbfontsz)
rc("ytick", labelsize=lbfontsz)
rc("axes", linewidth=lwidth)
#figA.clear()
matplotlib.pyplot.rcParams.update({'axes.titlesize': lbfontsz})
matplotlib.pyplot.rcParams.update({'font.size': lbfontsz})
matplotlib.pyplot.rcParams.update({'lines.linewidth': lwidth})
matplotlib.pyplot.rcParams.update({'ytick.major.width': lwidth})
matplotlib.pyplot.rcParams.update({'xtick.major.width': lwidth})
matplotlib.pyplot.rcParams.update({'ytick.major.size': 10.0})
matplotlib.pyplot.rcParams.update({'xtick.major.size': 10.0})
matplotlib.pyplot.rcParams.update({'font.family': 'serif'})
#matplotlib.pyplot.rcParams.update({'font.serif': 'Times'})


matplotlib.pyplot.subplots_adjust(wspace=0, hspace=0)

plotAc=figA.add_subplot(2,1,1,yscale='linear') 
plotAd=figA.add_subplot(2,1,2)      #
	

Omega_bar=2*pi*nu*sqrt(2*R_e**3/R_g)/c
print('_O_^2',Omega_bar**2,'_O_',Omega_bar)
flattening=(0.788-0.515*R_g/R_e)*Omega_bar**2 
print(R_e*(1-flattening))


def Beloborodov(cos_psi):
    """Beloborodov's approximation for cos_alpha(cos_psi) light bending function
    takes the cos psi 
    returns the cos alpha and its derivative
    """
    return 1. + (cos_psi - 1.)/redshift**2 ,1./redshift**2

def Poutanen(u,y):
    return ( 1 - u )*y*( 1 + u*u*y*y/112 - np.e/1e2*u*y*( np.log( 1 - y/2 ) + y/2 ) )

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
                dlag=dpsi*b/c/(1+sr) #*2/2
                psi+= dpsi
                lag+= dlag
    return psi,lag
# flattening=0


sin_i=sin(incl)
cos_i=cos(incl)


for p in range(NSpots):
    sin_theta=sin(theta[p])
    cos_theta=cos(theta[p])

    R=R_e*(1 - flattening*cos_theta**2) 
    dR=2*R_e*flattening*cos_theta*sin_theta # dR / d\theta

    u = R_g/R
    redshift=1.0/sqrt(1.0 - u) # 1/sqrt(1-R_g/R) = 1+ z = redshift
    f=redshift/R*dR
    sin_gamma=f/sqrt(1 + f**2) # angle gamma is positive towards the north pole 
    cos_gamma=1.0/sqrt(1 + f**2)
    beta=2*pi*nu*R*redshift*sin_theta/c
    Gamma=1.0/sqrt(1.0 - beta**2)
    Gamma1= (1.0-sqrt(1.0 - beta**2) )/ beta
    print('theta: ',theta[p],'gamma: ',arctan2(f,1.0))
     

    for t in range(NPhi):
          # if True: # find mu
                phi0 = phi[t]*pi/180
                sin_phi=sin(phi0)
                cos_phi=cos(phi0)
                cos_psi=cos_i*cos_theta + sin_i*sin_theta*cos_phi
                sin_psi=sqrt(1. - cos_psi**2)
               
                cos_alpha = 1.0 - Poutanen(u, 1.0 - cos_psi) # insert exact formula here
                sin_alpha = sqrt(1. - cos_alpha**2)
                sin_alpha_over_sin_psi= sin_alpha/sin_psi if sin_psi > 1e-4 else 1./redshift
                
                cos_xi = - sin_alpha_over_sin_psi*sin_i*sin_phi
                delta = 1./Gamma/(1.-beta*cos_xi)
                cos_sigma = cos_gamma*cos_alpha + sin_alpha_over_sin_psi*sin_gamma*(cos_i*sin_theta - sin_i*cos_theta*cos_phi)

                sin_sigma = sqrt(1. - cos_sigma**2)
                mu0=delta*cos_sigma # cos(sigma')
                # Omega=dS[p]*mu0*redshift**2*dcos_alpha #*Gamma*R*R/cos_gamma # 
         
                sin_chi_0= - sin_theta*sin_phi # times sin psi
                cos_chi_0=sin_i*cos_theta - sin_theta*cos_i*cos_phi # times sin psi 
                chi_0=arctan2(sin_chi_0,cos_chi_0)
                #  chi_0 = 0
                
                sin_chi_1=sin_gamma*sin_i*sin_phi*sin_alpha_over_sin_psi #times sin alpha sin sigma 
                cos_chi_1=cos_gamma - cos_alpha*cos_sigma  #times sin alpha sin sigma 
                chi_1=arctan2(sin_chi_1,cos_chi_1)
                
                sin_lambda=sin_theta*cos_gamma - sin_gamma*cos_theta
                cos_lambda=cos_theta*cos_gamma + sin_theta*sin_gamma
                cos_eps = sin_alpha_over_sin_psi*(cos_i*sin_lambda - sin_i*cos_lambda*cos_phi + cos_psi*sin_gamma) - cos_alpha*sin_gamma
                # this line is the longest one
                # alt_cos_eps=(cos_sigma*cos_gamma - cos_alpha)/sin_gamma # legit! thanks God I checked it!
                # sin_chi_prime=cos_eps*mu0*Gamma*beta # times something
                sin_chi_prime=cos_eps*mu0*delta*Gamma*beta*(1-Gamma1*cos_xi)# times something
                # cos_chi_prime=1. - cos_sigma**2 /(1. - beta*cos_xi) # times the samething

                cos_chi_prime=sin_sigma**2 - Gamma*mu0**2*beta*cos_xi*(1 - Gamma1*cos_xi)  # times the samething
                chi_prime=arctan2(sin_chi_prime,cos_chi_prime)   

                chi[t] = (chi_0+chi_prime+ chi_1)*180/pi#%(2*pi)

                cos_eps_sph = sin_alpha_over_sin_psi*(cos_i*sin_theta - sin_i*cos_theta*cos_phi)
                sin_chi_prime_sph=cos_eps_sph*mu0*delta*Gamma*beta*(1-Gamma1*cos_xi)# times something
                cos_chi_prime_sph=sin_alpha**2 - Gamma*mu0**2*beta*cos_xi*(1 - Gamma1*cos_xi)  # times the samething
                chi_prime_sph= arctan2(sin_chi_prime_sph,cos_chi_prime_sph)   

                chi_sph[t] = (chi_0+ chi_prime_sph)*180/pi#%(2*pi)



    plotAc.plot(phi,chi,"-",markersize=5,color=colors[p] )
    plotAc.plot(phi,chi_sph,"--",markersize=4,color=colors[p])
    plotAd.plot(phi,chi-chi_sph,"-",markersize=4,color=colors[p])


fontsize = 30
# plotAc.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
# plotAd.xaxis.set_major_formatter(matplotlib.pyplot.NullFormatter())
plotAd.set_xlabel(r'$\varphi\,[\degree]$',fontsize=fontsize)
plotAd.set_xticks([0,60,120,180,240,300,360])
plotAc.set_yticks([0,-90,-180,90,180])

plotAc.tick_params(axis='both', which='major', labelsize=fontsize)
plotAd.tick_params(axis='both', which='major', labelsize=fontsize)
# plotAd.set_ylabel(r'$\chi_{\mathrm{obl}}-\chi_{\mathrm{sph}}$',fontsize=fontsize)
plotAd.set_ylabel(r'$\chi_{{obl}}-\chi_{{sph}},\,[\degree]$',fontsize=fontsize)
plotAc.set_ylabel(r"$\chi_1 + \chi{'},\,[\degree]$",fontsize=fontsize)

# figA.tight_layout()

figA.subplots_adjust(left=0.15)

figA.savefig('fig2.pdf')#.format(e))
figA.clf()
