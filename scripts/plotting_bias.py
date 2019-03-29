import numpy as np
import matplotlib.pyplot as plt
import compute_rad_energy_from_waveforms as crefw
from lal import PI as LAL_PI
from lal import MTSUN_SI as LAL_MTSUN_SI
import plotsettings

date = '2015-02-08'

#M = [30., 50., 100.]
#Q = [1., 2., 4.]
#n_isco = [0.5, 0.75, 1.0]

M = [50.0]; Q = [2.0]; n_isco = [2.0, 2.5, 3.0, 3.5];


# Loading Data

delta = np.zeros(len(n_isco))
delta_Ef = np.zeros(len(n_isco))
sigma = np.zeros(len(n_isco))
snr = np.zeros(len(n_isco))


plt.figure(figsize=(16,8))
for i in range(len(M)):
  for j in range(len(Q)):
    for k in range(len(n_isco)):
#	M1 = M[i]/(1.+Q[j])
#	M2 = (M[i]*Q[j])/(1.+Q[j])

#	f_high = n_isco[k]*(1./6.)**(3./2)/(LAL_PI*M[i]*LAL_MTSUN_SI)

#	print f_high

	data_file = '../runs/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, str(M[i]), str(Q[j]), str(n_isco[k]))

	mtotal, logl, m1, m2, mc, q, eta, time = np.loadtxt(data_file, skiprows =1, usecols = (5, 7, 13, 16, 23, 29, 32, 34), unpack=True)

#	logl_max = max(logl)
#	snr[k] = np.sqrt(-2.*logl_max)
	
#	f_th, hf_th, ef_th, Ef_th = crefw.edf(M1, M2)
#	f_post, hf_post, ef_post, Ef_post = crefw.edf(np.mean(m1), np.mean(m2))


	delta[k] = abs(1041033614.-np.mean(time))/1041033614.
#	delta_Ef[k] = Ef_post[-1] - Ef_th[-1]
	sigma[k] = np.std(mtotal)/M[i]

    #plt.subplot(3, 3, 3*i+j+1)
    plt.plot(n_isco, delta, marker='o', color='r', label='gps time of injection used as trigger=1041033614.')
    #plt.plot(n_isco, sigma, marker='.', color='g', label='$\sigma$')
    plt.legend(loc='best')
    plt.title('$M=$'+str(M[i])+'; $q$='+str(Q[j]))
    plt.xlabel('$time$')
    plt.ylabel('$\Delta_{time}')    
    plt.grid(True)
    plt.hold(True)
plt.tight_layout()
"""
#	plt.subplot(1,3,1)
#	plt.plot(f_th, abs(hf_th), marker='o','r')
#	plt.plot(f_post, abs(hf_post), marker='o','g')
#	plt.subplot(1,3,2)
#	plt.plot(f_th, ef_th, marker='o','r')
#	plt.plot(f_post, ef_post, marker='o','g')
#	plt.subplot(1,3,3)
#	plt.plot(f_th, Ef_th, marker='o','r')
#	plt.plot(f_post, Ef_post,marker='o', 'g')


#f, hf, ef, Ef = crefw.edf(15., 15., 100.)

#plt.subplot(2,3,4)
#plt.plot(f, abs(hf), marker='o', 'r')
#plt.subplot(2,3,5)
#plt.plot(f, ef, marker='o','r')
#plt.subplot(2,3,6)
#plt.plot(f, Ef, marker='o','r')


plt.figure(figsize=(16,8))
#plt.plot(n_isco, delta_Ef)
plt.subplot(1,2,1)
plt.plot(n_isco, delta, marker='o', label='$\Delta$')
#plt.plot(n_isco, np.log10(snr),  marker='o', label='$log_{10}(SNR)$')
plt.plot(n_isco, 1./snr, marker='o', label='1/SNR')
plt.plot(n_isco, sigma,  marker='o', label='$\sigma$')
plt.legend(loc='best')
plt.xlabel('$f_{high}/f_{isco}$')
plt.ylabel('M (M_{\odot})')
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(n_isco, delta, marker='o', label='$\Delta$')
plt.plot(n_isco, 1./snr, marker='o', label='1/SNR')
plt.plot(n_isco, sigma,  marker='o', label='$\sigma$')
plt.plot(n_isco, sigma[3]/((snr/snr[3])), marker='o', label='1/SNR')
#plt.plot(n_isco, 1./(9.204*snr), marker='o', label='scaled 1/SNR')
#plt.loglog(1./snr, sigma, 'o-')
plt.legend(loc='best')
plt.xlabel('$f_{high}/f_{isco}$')
plt.ylabel('M (M_{\odot})')
plt.grid(True)

plt.suptitle('$M=$'+str(M[0])+'; $q$='+str(Q[0]))
"""
plt.show()
