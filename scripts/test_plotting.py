import numpy as np
import matplotlib.pyplot as plt
import compute_rad_energy_from_waveforms as crefw
from lal import PI as LAL_PI
from lal import MTSUN_SI as LAL_MTSUN_SI
import plotsettings
import commands
import my_plotter as mp
import loadline as ll
import os

date = '2015-02-25'
approximant = 'IMRPhenomB'

#phase_list = ['inspiral']
phase_list = ['ringdown']
M_list = [100.]
Q_list = [4.]
#n_isco_list = [1.0]
n_isco_list = [0.25]

# Loading Data

for phase in phase_list:
  for M in M_list:
    for Q in Q_list:
      for n_isco in n_isco_list:
        M1 = M/(1.+Q)
        M2 = M*Q/(1.+Q)
        ETA = (M1*M2)/((M1+M2)**2.)

        MF = M*(1. + (np.sqrt(8./9.)-1.)*ETA - 0.4333*(ETA**2.) - 0.4392*(ETA**3.))
        AF_MF = ETA*np.sqrt(12.) - 3.871*(ETA**2.) + 4.028*(ETA**3.)


        data_file = '../runs/%s/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, phase, str(M), str(Q), str(n_isco))
        legend = ll.load_line(data_file, dtype='string')

        mtotal, logl, m1, m2, mc, q, eta = np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('mtotal'), list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2'), list(legend).index('mc'), list(legend).index('q'), list(legend).index('eta')), unpack=True)

        mf_M = 1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.)
        af_mf = eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.)


	#H_af_mf, af_mf_edges = np.histogram(af_mf, bins=100)
	#H_mf_M, mf_M_edges = np.histogram(mf_M, bins=100)
	
	H, af_mf_edges, mf_M_edges = np.histogram2d(af_mf, mf_M, bins=100)

	f = open('data_%s_2d.dat'%phase, 'w')
	for i in range(len(af_mf_edges[:-1])):
	  for j in range(len(mf_M_edges[:-1])):
		f.write('%f %f %f \n'%(af_mf_edges[i], mf_M_edges[j], H[i,j])) 
	f.close()
	
	#np.savetxt('data_%s.dat'%phase, (mf_M_edges[:-1], H_mf_M, af_mf_edges[:-1], H_af_mf ))
	#plt.figure()
	#plt.hist(mf_M, bins=100)
	#plt.show()
	
	"""
	plt.figure()
	plt.subplot(2,2,1)
	plt.hist(af_mf, bins=75, color='Aquamarine', normed=True)
	plt.axvline(x=AF_MF, ls='--', color='g')
	plt.xlabel('$m_1$')
	plt.ylabel('$P(m_1)$')
	plt.xlim([min(af_mf), max(af_mf)])
	plt.grid()	
	
	plt.subplot(2,2,3)
	plt.hist2d(af_mf,mf,bins=75, cmin= 0.00001, cmap='hot', normed=True)
	plt.axhline(y=MF, lw=2., ls='--', color='r')
	plt.axvline(x=AF_MF, lw=2., ls='--', color='r')
	plt.xlim([min(af_mf), max(af_mf)])
	plt.ylim([min(mf), max(mf)])
	plt.ylabel('$m_2$')
        plt.xlabel('$m_1$')
	plt.grid()

	plt.subplot(2,2,4)
	plt.hist(mf, bins=75, normed=True, color='Aquamarine', orientation='horizontal')
	plt.axhline(y=MF, ls='--', color='g')   
        plt.ylabel('$m_2$')
        plt.xlabel('$P(m_2)$')
	plt.ylim([min(mf), max(mf)])
        plt.grid()
	plt.tight_layout()
plt.show()
"""
