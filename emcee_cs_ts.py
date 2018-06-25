#Running emcee for polarization angles

import numpy as np
from numpy import logspace, zeros, fromfile
from numpy import pi, exp, log, sqrt, sin, cos, arccos, arctan2
#print (np.__file__)
#import numpy.random as random
#import scipy
try:
	import matplotlib
	matplotlib.use('Agg')
	##import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt
	import corner
	import kombine
	from pyf_read_data import pyf_read_data
except ImportError:
	pass
import emcee
#from emcee import PTSampler
#from emcee.utils import MPIPool
#from emcee import EnsembleSampler
import math, os
#import pymultinest
import time
from cs_ts2_func import compf

spath = "res/B/"#str(sys.argv[1])

NEnergy = 281
NPhase = 150
x_l, x_u = -3.7 , .3 # lower and upper bounds of the log_10 energy span
evere=.5109989e6 # electron volts in elecron rest energy 
IntEnergy = logspace(x_l,x_u,NEnergy), log(1e1)*(x_u-x_l)/(NEnergy-1.) 
x,x_weight=IntEnergy #energies
ene_index = 118
energy_keV = x*evere/1e3 #energies in keV
#print(energy_keV[110:160])
print("The energy chosen = ", energy_keV[ene_index], " keV")
param_names = ["mass","rad","incl","theta"]
params_true = [1.4,12.0,40.0,60.0]
low_limit = [1.0, 4.0, 0.0, 0.0]
high_limit = [2.0, 18.0, 90.0, 90.0]
stlow_limit = np.copy(low_limit)
sthigh_limit = np.copy(high_limit)
restart = False


def readdata():
	PulsName='res/B/B0P2'
	inFlux = open(PulsName+'FF.bin')
	inphi = open(PulsName+'ff.bin')
	Flux1 = fromfile(inFlux)
	phi = fromfile(inphi)
	fluxlcurve_Iene = Flux1[0+ene_index*3:len(Flux1):3*NEnergy]
	fluxlcurve_Qene = Flux1[1+ene_index*3:len(Flux1):3*NEnergy]
	fluxlcurve_Uene = Flux1[2+ene_index*3:len(Flux1):3*NEnergy]
	Flux2 = np.array([fluxlcurve_Iene, fluxlcurve_Qene, fluxlcurve_Uene])
	return phi, Flux2


#likelyhoods for emcee:
def lnprob(modelpar, low_limit, high_limit):

	#check boundaries:
	for ii in range(0,len(modelpar)):
		if(modelpar[ii]<low_limit[ii] or modelpar[ii] > high_limit[ii]):
			#print(modelpar[ii],low_limit[ii],high_limit[ii])
			#quit()
			return -np.inf

	mass = modelpar[0]
	rad = modelpar[1]
	incl = modelpar[2]
	theta = modelpar[3]

	if(mass/rad > 0.9*1.0/(2.95*1.52)): #checking causality                
		return -np.inf

	Flux = compf(mass,rad,incl,theta)
	#print(Flux)
	phi,Flux_obs = readdata()

	ene = ene_index #the chosen energy
	I=zeros(NPhase)
	Q=zeros(NPhase)
	U=zeros(NPhase)
	for t in range(NPhase):
		I[t],Q[t],U[t]=Flux[t,ene]
	#p=sqrt(Q**2+U**2)/I*100
	PA=arctan2(-U,-Q)*90/pi+90
	for t in range(NPhase):
		I[t],Q[t],U[t]=Flux_obs[0,t], Flux_obs[1,t], Flux_obs[2,t]
	PA_obs=arctan2(-U,-Q)*90/pi+90
	
	#print(PA)
	#print(PA_obs)
	#quit()

	#sigma_tot2 = modchan[idat][iene]*expo[idat]+pow(insigma,2)+std::pow(modchan[idat][iene]*expo[idat]*0.005,2); //pearson (model error) 
	#//norm = 0.5*log(2.0*pi*sigma_tot2); //when comparing to normalized poisson
	#norm = 0.5*log(sigma_tot2);
	#if(cnts[idat][iene] > 20){
	#    logprobsum = logprobsum-smbeta*pow((cnts[idat][iene]-modchan[idat][iene]*expo[idat]),2)/(2.0*sigma_tot2)-norm;
	#}
	loglik = 0.0
	insigma = 0.0
	for t in range(NPhase):
		sigma_tot2 = 1.0#PA+insigma**2+(0.005*PA)**2
		norm = 0.0#0.5*log(sigma_tot2)
		loglik = loglik - (PA[t]-PA_obs[t])**2/(2.0*sigma_tot2)-norm
	#print("loglikelihood = ",loglik)
	#print(PA[0],PA_obs[0])
	return loglik


#Flux = compf(1.0,1.0,1.0,1.0)
print("Running emcee")

nwalkers = 20#14#50
nparams = len(param_names)
ndim = nparams
p0 = np.zeros((nwalkers,ndim))


#test call to lnprob:
start = time.time()
print(lnprob(params_true,low_limit,high_limit))
end = time.time()
print ("Time spend for one fit: ")
print(end - start)
#quit()

jj = 0
# Choose an initial set of positions for the walkers.
#p0 = [np.random.rand(ndim) for i in xrange(nwalkers)]
#print(p0[:][0])
while(jj < nwalkers):
	#print jj
	causal = False
	for ii in range(0,ndim):
		p0[jj][ii] = stlow_limit[ii]+np.random.rand()*(sthigh_limit[ii]-stlow_limit[ii])
	#checking additional constraints: (if not fullfiled, guess new starting point to all params)
	if(p0[jj][0]/p0[jj][1] > 0.96*1.0/(2.95*1.52)): #checking causality
		continue
	jj = jj+1


print("walkers initialized succesfully")
#print(p0)

## Initialize the sampler with the chosen specs.
start = time.time()
print("Initializing the sampler")
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[low_limit,high_limit])#,threads=nthreads, pool=pool)#, args=[means, icov])
#note: threads ignored when pool used for parallelization
#set pool to None if want to use multiprocessing parallelization
end = time.time()
print ("Time spend: ")
print(end - start)

f = open(spath + "emcee.dat", "w")
f.close()
start = time.time()
## Run some steps as a burn-in.
#pos, prob, state = sampler.run_mcmc(p0, 50)
nsteps = 50#2#3#2#50

#result = sampler.sample(p0, iterations=nsteps)

#sampler.run_mcmc(p0, nsteps)

#print("hello??")
#quit()

print("Start burn-in computation:")
for i, result in enumerate(sampler.sample(p0, iterations=nsteps)): #or use sampler.run()?:
	print(i)	
	if (i+1) % 10 == 0:
		print("{0:5.1%}".format(float(i) / nsteps))
		print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
	pos, prob, state = result
	position = pos #result[0]
	#print(position)
	f = open(spath + "emcee.dat", "a")
	for k in range(position.shape[0]):
		f.write("{0:4d} {1:s}".format(k, " ".join(map(str,position[k]))))
		f.write(" ")
		f.write("{:f}\n".format(prob[k]))
	f.close()
	if(i == 0):
		end = time.time()
		print ("Time spend on the first step: ")
		print(end - start)


## Reset the chain to remove the burn-in samples.
print("Reset chain to remove the burn-in samples")
sampler.reset()


## Starting from the final position in the burn-in chain, sample for 100
## steps.
#sampler.run_mcmc(pos, 1000, rstate0=state) #rstate0 = The state of the random number generator.

nsteps = 4000#2#20000
for i, result in enumerate(sampler.sample(pos, iterations=nsteps, rstate0=state)): #or use sampler.run()?:
	if (i+1) % 100 == 0:
		print("{0:5.1%}".format(float(i) / nsteps))
		print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
	pos, prob, state = result
	position = pos #result[0]
	#print(position)
	f = open(spath + "emcee.dat", "a")
	for k in range(position.shape[0]):
		f.write("{0:4d} {1:s}".format(k, " ".join(map(str,position[k]))))
		f.write(" ")
		f.write("{:f}\n".format(prob[k]))
	f.close()


##samples = sampler.chain[:, 50:, :].reshape((-1,ndim))
samples = sampler.chain[:, 0:, :].reshape((-1,ndim))

print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

limits =  zip(low_limit,high_limit)
fig = corner.corner(samples,labels=param_names,truths=params_true,range=limits)
fig.savefig(spath+"emcmc_triangle.pdf")

#print lnprob()









