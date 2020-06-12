###############################
#This program can be used to calculate polarized pulse profiles from accreting millisecond pulsars
###############################

<<<<<<< HEAD
#from cs_ts2_comp_pa import compf

mass = 1.4 # 1.4
rad = 12.0 # 12.0
incl = 40.0 #60.0 #60.0 #40.0 # 40.0
theta = 60.0#20.0 #20.0 #60.0 #-120.0#60.0 # 60.0
rho = 10.0 #1.0 #1.0 # 10.0
freq = 600.0

#mass and rad for non-rotating star if req=12 and M=1.4 for the rotating one (see mr_rot_nonrot.py):
#mass, rad = 1.39404589, 11.42791781

#Flux = compf(mass,rad,incl,theta,rho,freq,spherical=False)

from cs_ts2_func import compf
Flux = compf(mass,rad,incl,theta,rho,spherical=True)


loop = False#True
if(loop):

	incls = [50.0]#[40.0,60.0,80.0]
	thetas = [90.0]#[10.0,30.0,60.0,90.0]#[10.0,20.0,40.0,60.0]
	freqs = [600.0]#[1.0,100.0,300.0,600.0]
	sph = [True]#[False,True]
=======
from cs_ts2_comp_pa import compf

mass = 1.4 # 1.4
rad = 12.0 # 12.0
incl = 60.0 #60.0 #40.0 # 40.0
theta = 20.0 #20.0 #60.0 #-120.0#60.0 # 60.0
rho = 1.0 #1.0 # 10.0
freq = 600.0

#Flux = compf(mass,rad,incl,theta,rho,freq,spherical=False)
loop = True
if(loop):

	incls = [40.0,60.0,80.0]
	thetas = [10.0,20.0,40.0,60.0]
	freqs = [1.0,100.0,300.0,600.0]
	sph = [False,True]
>>>>>>> 8396ea969cf96e78354ff847d3eae8b8bbf87d06
	for ifr in range(len(freqs)):
		for i in range(len(incls)):
			for it in range(len(thetas)):
				for isp in range(len(sph)):
					freq = freqs[ifr]
					incl = incls[i]
					theta = thetas[it]
					spherb = sph[isp]
					#spath = "pulses/pd_yes0_b_no0/pulse_thom"+str(int(antpd)+1)+"_r"+str(int(rho))+"t"+str(int(theta))+"i"+str(int(incl))+"p"+str(int(10000*pol))
					#if(incl!=60.0 and theta!=20.0):
					#	continue
					#Flux = compf(mass,rad,incl,theta,rho,pol,spherical=False,antipodal=antpd,spath=spath)
					Flux = compf(mass,rad,incl,theta,rho,freq,spherical=spherb)
	exit()


