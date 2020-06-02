###############################
#This program can be used to calculate polarized pulse profiles from accreting millisecond pulsars
###############################

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


