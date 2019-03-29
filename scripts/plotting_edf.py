import numpy as np
import matplotlib.pyplot as plt
import compute_rad_energy_from_waveforms as crefw
import os
from lal import PI as LAL_PI
from lal import MTSUN_SI as LAL_MTSUN_SI
import plotsettings
import commands
import my_plotter as mp
import loadline as ll


date = '2015-02-16'

f_low = 10.
fs = 2048.
f_final = fs/2.

M = [100.] #[30., 50., 100.]
Q = [1., 2., 4.]
n_isco = [0.25, 0.33, 0.5, 0.8, 0.9, 1.0]

delta_Ef = np.zeros(len(n_isco))
sigma = np.zeros(len(n_isco))


plt.figure(figsize=(16,4))
for i in range(len(M)):
  for j in range(len(Q)):
    for k in range(len(n_isco)):
	f_high = n_isco[k]*(1./6.)**(3./2)/(LAL_PI*M[i]*LAL_MTSUN_SI)

	data_file = '../runs/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, str(M[i]), str(Q[j]), str(n_isco[k]))
	legend = ll.load_line(data_file, dtype='string')

	mtotal, logl, m1, m2, mc, q, eta = np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('mtotal'), list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2'), list(legend).index('mc'), list(legend).index('q'), list(legend).index('eta')), unpack=True)


	M1 = M[i]/(1.+Q[j])
        M2 = M[i]*Q[j]/(1.+Q[j])
	

	f_th, hf_th, ef_th, Ef_th = crefw.edf(M1, M2, f_low, fs, f_final)
	f_post, hf_post, ef_post, Ef_post = crefw.edf(np.mean(m1), np.mean(m2), f_low, fs, f_final)
	
	for l in range(len(f_post)):
		if f_post[l] > f_high:
                        Ef_wt_post = Ef_post[l]
			break

	delta_Ef[k] = Ef_wt_post/Ef_th[-1]
	sigma[k] = np.std(mtotal)/M[i]

    plt.subplot(1,3,j+1)
    plt.loglog(n_isco, delta_Ef, marker='o', label='$E(f_{wt})/E_{tot, inj}$')
    plt.loglog(n_isco, sigma, marker='o', label='$\sigma _{M_{total}/M_{inj}}$')
    plt.legend(loc='best')
    plt.grid(True, which='both')
    plt.xlabel('$n_{isco}$')
    plt.title('M='+str(M[i])+'; Q='+str(Q[j]))
    plt.tight_layout()
plt.savefig('../runs/%s/f_high_plot_%s.png'%(date, str(M[i])))   
#plt.show()
