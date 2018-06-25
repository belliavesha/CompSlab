#import:
import matplotlib
matplotlib.use('Agg')
from numpy import linspace, logspace, empty, zeros, ones, array, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
from numpy import absolute, sign, floor, ceil
from numpy.polynomial.laguerre import laggauss
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import interp1d
from scipy.special import kn
from matplotlib.pyplot import *
from bisect import bisect


#physical constants:
evere=.5109989e6 # electron volts in elecron rest energy 
G=13275412528e1 # G*M_sol in km^3/s^2 
c=299792458e-3 # speed of light in km/s

NPhase = 150#128*2#120#128#120 # Number of equidistant phase points
NBend= 20 # Number of knots in light bending integrations
NBendPhase= 1000 # Number of psi/aplha grid points
IntBend = leggauss(NBend)
NEnergy = 281 # 50# 101 # number of energy points (x)

phi,phi_weight=linspace(0,2*pi,num=NPhase,endpoint=False,retstep=True) #Size of spacing between samples = phi_weight
nu=600 # star rotation frequency in Hz
M=1.4 # star mass in solar masses
R_g=M*2.95 # gravitational Schwarzschild radius
R_e=12.0 # equatorial radius of the star in kilometers      
i=pi*4/18#1.5    # line of sight colatitude #in radians
#i=pi/180.0*80.0
theta = [pi/3,pi-pi/3] # spot colatitude
#theta = [pi/180.0*15.0,pi-pi/180.0*15.0] # spot colatitude

#Try switching:
#theta_new = [i,pi-i] # spot colatitude
#i = theta[0] # spot colatitude
#theta = theta_new


tau_T= .6 # Thomson optical depth of thermalization 
x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
Theta = 0.1  # dimensionless electron gas temperature (Theta = k T_e / m_e c^2) # it's about 0.1 
T = 0.002 # 10/evere #  dimensionless photon black body temperature T = k T_bb / m_e c^2
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) # sample points and weights for integrations over the spectrum computing sorce function
x,x_weight=IntEnergy


shapes = ["Sphere","AlGendy"]

colors = ["green","blue"]

figA=figure(2,figsize=(10,14))
figA.clear()

plot_only_I = False#True
if(plot_only_I):
	figA=figure(1,figsize=(14,14))
	figA.clear()
	plotAF=figA.add_subplot(1,1,1,yscale='linear') 
else:	
	plotAF=figA.add_subplot(3,1,1,yscale='linear') 
	plotAp=figA.add_subplot(3,1,2)      #
	plotAc=figA.add_subplot(3,1,3)      #

for ish in range(1,len(shapes)):

	#oblateness='AlGendy'#'Sphere'#'AlGendy'
	oblateness=shapes[ish]

	#AtmName='res/B/C1obl' # the prefix for all result files related to the set of parameters
	#AtmName='res/B/B0P2' # the prefix for all result files related to the set of parameters
	if(ish == 0):
		#AtmName='res/C2/C1sph'
		#PulsName='res/C2/sph2test'
		PulsName='res/B/B0P2'
	else:
		#AtmName='res/C2/C1obl'  
		#PulsName='res/C2/obl2test'
		PulsName='res/B/B0P2'
	#PulsName=AtmName+'P1'
	computePulse= True
	plotAtm=not True
	plotPulse=True
	mod=True

	#outF = open(PulsName + 'F.bin','w')
	#outf = open(PulsName + 'f.bin','w')
	#Flux.tofile(outF,format="%e")
	#phi.tofile(outf,format="%e")


	inFlux = open(PulsName+'FF.bin')
	inphi = open(PulsName+'ff.bin')
	Flux1 = fromfile(inFlux)
	phi = fromfile(inphi)
	#print(phi, Flux1)	
	fluxlcurve0 = Flux1[0:len(Flux1):3*NEnergy] #light curve with lowest E
	fluxspec0 = Flux1[0:3*NEnergy:3] #spectrum at phase=0
	#print(fluxlcurve0)
	#print(" ")
	#print(fluxspec0)

	ene = 166#140#166#140 #The chosen energy index
	fluxlcurve_Iene = Flux1[0+ene*3:len(Flux1):3*NEnergy]
	fluxlcurve_Qene = Flux1[1+ene*3:len(Flux1):3*NEnergy]
	fluxlcurve_Uene = Flux1[2+ene*3:len(Flux1):3*NEnergy]
            	            
	#Flux=zeros((NPhase,NEnergy,3))
	#print(fluxlcurve_Iene)
	#print(fluxlcurve_Qene)
	#print(fluxlcurve_Uene)
	
	labelsize=20
	fontsize=25

	#phase=list(phi/2/pi)+[1.]
	phase=list(phi)+[1.]
	I=zeros(NPhase+1)
	Q=zeros(NPhase+1)
	U=zeros(NPhase+1)
	for t in range(NPhase+1):
		#I[t],Q[t],U[t]=Flux[t-1,e]*x[e] 
		I[t],Q[t],U[t]=fluxlcurve_Iene[t-1]*x[ene] ,fluxlcurve_Qene[t-1]*x[ene] ,fluxlcurve_Uene[t-1]*x[ene] 

	p=sqrt(Q**2+U**2)/I*100
	PA=arctan2(-U,-Q)*90/pi+90
	#PA=arctan2(U,Q)*180/pi#+90

	#if(ish > 0):
	#	plotAF.plot(phase,I/I.max())#,color=col)
	#	plotAp.plot(phase,p)#,color=col)
	#	plotAc.plot(phase,PA)#,color=col)	 
	#	break	           	

	figA.suptitle(r'$\nu={:5.0f}Hz$'.format(nu)+
	              r'$,\,R_e={:5.1f}km$'.format(R_e)+
	              r'$,\,M=$'+str(M)+r'$M_{\odot}$'+',\n'+
	              r'$\,\theta={:5.1f}\degree$'.format(theta[0]*180/pi)+
	              r'$,\,i={:5.1f}\degree$'.format(i*180/pi)+',\n'+
	              r'$\,x[keV]={:6.2f}$'.format(x[ene]*evere/1e3),fontsize=fontsize)  

	# plotAF.set_xlim([x[0],x[-1]])
	plotAF.set_xlim(0,1)

	# plotAF.locator_params(axis='y', nbins=10)
	plotAF.set_ylabel(r"$F_{x}(\varphi)/F_{x}^{\mathrm{max}}$",fontsize=fontsize)
	plotAF.tick_params(axis='both', which='major', labelsize=labelsize)
	# plotAF.plot(x,xIinx,'k-.')

	# plotAF.set_xlim([x[0],x[-1]])
	if not(plot_only_I):
		plotAp.set_xlim(0,1)
		plotAp.tick_params(axis='both', which='major', labelsize=labelsize)
		plotAp.set_ylabel(r'$p\,[ \% ]$',fontsize=fontsize)
		plotAc.set_xlim(0,1)
		plotAc.set_ylim(0,180)
		plotAc.set_yticks([0,30,60,90,120,150,180])
		#plotAc.set_ylim(-180,180)
		#plotAc.set_yticks([0,30,60,90,120,150,180])
		plotAc.tick_params(axis='both', which='major', labelsize=labelsize)
		plotAc.set_ylabel(r'$\chi\,[\degree]$',fontsize=fontsize)
		plotAc.set_xlabel(r'$\varphi\,[360\degree]$',fontsize=fontsize)

	#col=colors[(e*NColors)//NEnergy]
	col = colors[ish]

	plotAF.plot(phase,I/I.max(),color=col)
	#plotAF.plot(phase,I,color=col)
	if not(plot_only_I):
		plotAp.plot(phase,p,color=col)
		plotAc.plot(phase,PA,color=col)

compare_to_tslc = False#True
if(compare_to_tslc):
	colors = ["orange","red"]
	for ic in range(1,2):

		datafile=""
		if(ic == 0):
			datafile = "../../MCMC/results_sph_test.txt"
		if(ic == 1):
			datafile = "../../MCMC/results_obl2_test.txt"
		input = file(datafile, 'r')
		lines = input.readlines()
		input.close()

		Nchain_size = 256#128#16
		#c_lines = 21#28
		c_lines = 29
		egrid = 40
		full_chain= [[] for x in xrange(egrid+1)]

		input = file(datafile, 'r')
		lines = input.readlines()
		input.close()

		for j in range(0,len(full_chain)):
			for i in range(c_lines,c_lines+Nchain_size): #not reading comment lines
				parts = lines[i].split()
				full_chain[j].append(float(parts[j]))
			parts = lines[c_lines].split()

		full_chain = np.array(full_chain)

		#print full_chain[0,:]
		#print full_chain[30,:]
		energy_keV = np.arange(2.0,25.00001, (25.0-2.0)/int(egrid-1))
		phase = full_chain[0,:]
		norm_obsF = np.zeros((len(phase), egrid))
		for i in range(1,egrid+1):
			norm_obsF[:,i-1] = full_chain[i,:]*energy_keV[i-1]/np.max(full_chain[i,:]*energy_keV[i-1])


		#for i in range(0,egrid):
		ene = 14#38##14
		print("energy_keV = ", energy_keV[ene])

		phshift = 0.0#-0.005#0.025
		for i in range(ene,ene+1):
			#frame1 = plt.gca()
			#line = frame1.axes.plot(phase,norm_obsF[:,i],markersize=2,label=str("%.2f" % energy_keV[i])+' keV',color=colors[ic])
			plotAF.plot(phase-phshift,norm_obsF[:,i],markersize=2,label=str("%.2f" % energy_keV[i])+' keV',color=colors[ic])




#figA.savefig('res/C2/obl_sph_comp.pdf')#.format(e))
figA.savefig('res/B/plot.pdf')#.format(e))
figA.clf()



