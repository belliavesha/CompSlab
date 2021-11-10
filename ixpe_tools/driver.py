###############################
#This program can be used to calculate and save polarized pulse profiles from accreting millisecond pulsars
###############################
import time
from polpulse import compf

mass = 1.4 # 1.4
rad = 12.0 # 12.0
incl = 60.0 #60.0 #40.0 # 40.0
theta = 20.0 #20.0 #60.0 #-120.0#60.0 # 60.0
rho = 1.0 #1.0 # 10.0
pol = 0.0 #0.1171 #0.0 #0.1171 #0.0 #0.1171 #0.0 #0.1171 # 0.1171 #USE NOW ONLY 0.0 or 0.1171


from numpy import logspace, zeros, fromfile, linspace
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
NEnergy = 50 
NPhase = 150
x_l, x_u = -3.7 , -1.2 # lower and upper bounds of the log_10 energy span 
evere=.5109989e6 # electron volts in elecron rest energy 
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) 
x,x_weight=IntEnergy #energies
phase =linspace(0,1,num=NPhase,endpoint=True,retstep=False) #input phase points
energy_keV = x*evere/1e3 # input energies in keV

loop = True
#if len(sys.argv) > 1:
if(loop):
	#print("Choosing param values for i,theta,rho from CL arguments:")
	#incl = float(sys.argv[1])
	#theta = float(sys.argv[2])
	#rho = float(sys.argv[3])
	#antpd = bool(sys.argv[4])
	#spath = str(sys.argv[5])

	print("Running in loop mode! Make sure not overwriting old results.")        
	incls = [60.0]#[50.0,60.0,70.0]
	thetas = [20.0]# [120.0]#[80.0, 100.0, 120.0,140.0,160.0]#[20.0]#[10.0,20.0,30.0]
	rhos = [1.0]#[30.0]#[1.0]
	antipds = [False]#[True,False]
	for ir in range(len(rhos)):
		for i in range(len(incls)):
			for it in range(len(thetas)):
				for ai in range(len(antipds)):
					rho = rhos[ir]
					incl = incls[i]
					theta = thetas[it]
					antpd = antipds[ai]
					spath = "pulses/pulse_thom"+str(int(antpd)+1)+"_r"+str(int(rho))+"t"+str(int(theta))+"i"+str(int(incl))+'_test_E0_new'
					#spath = "pulses/test/pulse_test_thom"+str(int(antpd)+1)+"_r"+str(int(rho))+"t"+str(int(theta))+"i"+str(int(incl))+"p"+str(int(10000*pol))
					#if(incl!=60.0 and theta!=20.0):
					#	continue
					start = time.time()
					Flux = compf(mass,rad,incl,theta,rho,pol,energy_keV,phase,spherical=False,antipodal=antpd,spath=spath,savePulse=True,atmos_path="atmos_thom/")
					end = time.time()
					print("Time pulse computation:",end - start)
					
					
else:
	Flux = compf(mass,rad,incl,theta,rho,pol,spherical=False,antipodal=False,spath="pulses/pulse_testX",savePulse=True,atmos_path="atmos_thom/") 

#import subprocess
#subprocess.call("./pulse_rename.sh",shell=True)
