import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import compute_rad_energy_from_waveforms as crefw
from lal import PI as LAL_PI
from lal import MTSUN_SI as LAL_MTSUN_SI
import plotsettings
import commands
import my_plotter as mp
import loadline as ll
import matplotlib.tri as mtri

date = '2015-02-21'

M = [100.] #[30., 50., 100.]
Q = [4.] #[1., 2., 4.]
n_isco = [0.7, 1.0]

color = np.genfromtxt('/archive/home/abhirup/Documents/color.dat', dtype='str')

# Loading Data

plt.figure(figsize=(16,len(n_isco)*2))
for i in range(len(M)):
  for j in range(len(Q)):
    for k in range(len(n_isco)):
	M1 = M[i]/(1.+Q[j])
        M2 = M[i]*Q[j]/(1.+Q[j])
	ETA = (M1*M2)/((M1+M2)**2.)
	
	MF = M[i]*(1. + (np.sqrt(8./9.)-1.)*ETA - 0.4333*(ETA**2.) - 0.4392*(ETA**3.))
        AF_MF = ETA*np.sqrt(12.) - 3.871*(ETA**2.) + 4.028*(ETA**3.)

	
	#data_file = '../runs/%s/inspiral/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, str(M[i]), str(Q[j]), str(n_isco[k]))
	data_file = '../runs/%s/ringdown/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, str(M[i]), str(Q[j]), str(n_isco[k]))
	legend = ll.load_line(data_file, dtype='string')

	mtotal, logl, m1, m2, mc, q, eta = np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('mtotal'), list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2'), list(legend).index('mc'), list(legend).index('q'), list(legend).index('eta')), unpack=True)

	mf = mtotal*(1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.))
	af_mf = eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.)

	plt.subplot(1,len(n_isco),k+1)
	plt.scatter(eta, mtotal, c=logl, vmax = max(logl), vmin = max(logl)-1, edgecolor='')
	plt.axvline(x=ETA, ls='--', lw=0.5)
        plt.axhline(y=M[i], ls='--', lw=0.5)

	plt.ylabel('$M_{total}$')
	plt.xlabel('$\eta$')        
	plt.colorbar()
	plt.title('n='+str(n_isco[k]))

plt.tight_layout()
#plt.savefig('../plots/%s_(SEOBNRv2)/inspiral/likelihood_%s_%s.png'%(date, str(M[i]), str(Q[j])))
plt.savefig('../plots/%s_(SEOBNRv2)/ringdown/likelihood_%s_%s.png'%(date, str(M[i]), str(Q[j])))
plt.show()
