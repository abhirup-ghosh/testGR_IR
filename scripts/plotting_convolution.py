import numpy as np
import scipy as sp
from scipy import signal
import matplotlib.pyplot as plt
import loadline as ll


date = '2015-02-25'
approximant = 'IMRPhenomB'
phase = ['inspiral', 'ringdown']
M, Q, n_isco, n_qnm = 100., 4., 1.0, 0.25

data_file = '../runs/%s/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, phase[0], str(M), str(Q), str(n_isco))
legend = ll.load_line(data_file, dtype='string')

mtotal, logl, m1, m2, mc, q, eta = np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('mtotal'), list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2'), list(legend).index('mc'), list(legend).index('q'), list(legend).index('eta')), unpack=True)

mf_i = 1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.)
af_mf_i = eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.)

data_file = '../runs/%s/%s/%s_%s/nisco_%s/cbc_bayes/posterior_samples.dat'%(date, phase[1], str(M), str(Q), str(n_qnm))
legend = ll.load_line(data_file, dtype='string')

mtotal, logl, m1, m2, mc, q, eta = np.loadtxt(data_file, skiprows =1, usecols = (list(legend).index('mtotal'), list(legend).index('logl'), list(legend).index('m1'), list(legend).index('m2'), list(legend).index('mc'), list(legend).index('q'), list(legend).index('eta')), unpack=True)

mf_r = 1. + (np.sqrt(8./9.)-1.)*eta - 0.4333*(eta**2.) - 0.4392*(eta**3.)
af_mf_r = eta*np.sqrt(12.) - 3.871*(eta**2.) + 4.028*(eta**3.)


#1D CONVOLUTIONS

H_mf_i, mf_edges_i = np.histogram(mf_i, bins=100)
H_mf_r, mf_edges_r = np.histogram(mf_r, bins=100)

H_af_i, af_edges_i = np.histogram(af_mf_i, bins=100)
H_af_r, af_edges_r = np.histogram(af_mf_r, bins=100)



H_mf = sp.interpolate.UnivariateSpline(mf_edges_i[:-1],H_mf_i)#, bbox=[mf_edges_i[0], mf_edges_i[-2]])
H_mf_i_int = H_mf(mf_edges_r[:-1])

H_af = sp.interpolate.UnivariateSpline(af_edges_i[:-1],H_af_i)#, bbox=[mf_edges_i[0], mf_edges_i[-2]])
H_af_i_int = H_af(af_edges_r[:-1])

plt.figure(figsize=(16,8))
#plt.plot(mf_edges_i[:-1], H_mf_i, ms='.', color='g')
#plt.plot(mf_edges_r[:-1], H_mf_r, color='r')
#plt.plot(mf_edges_r[:-1], H_mf_i_int, color='b')


delta_mf = 0.001
c_mf = np.convolve(H_mf_i_int, H_mf_r)*delta_mf

delta_af = 0.001
c_af = np.convolve(H_af_i_int, H_af_r)*delta_af


plt.subplot(2,2,4)
plt.plot(c_mf, np.linspace(-len(c_mf)/2, len(c_mf)/2, len(c_mf))*delta_mf, marker='o', ms=0.5)
plt.ylabel('$\Delta M_f$')
plt.xlabel('$P(\Delta M_f)$')
plt.grid()
plt.subplot(2,2,1)
plt.plot(np.linspace(-len(c_af)/2, len(c_af)/2, len(c_af))*delta_af, c_af, marker='o', ms=0.5)
plt.xlabel('$\Delta a_f/M_f$')
plt.ylabel('$P(\Delta a_f/M_f)$')
plt.grid()


#2D CONVOLUTIONS

H_i, af_edges_i, mf_edges_i = np.histogram2d(af_mf_i, mf_i, bins=50)
H_r, af_edges_r, mf_edges_r = np.histogram2d(af_mf_r, mf_r, bins=50)

H = sp.interpolate.RectBivariateSpline(af_edges_i[:-1], mf_edges_i[:-1], H_i)
H_int_i = H(af_edges_r[:-1], mf_edges_r[:-1])

c = sp.signal.convolve2d(H_int_i, H_r)*delta_af*delta_mf

af_len, mf_len = np.shape(c)[0], np.shape(c)[1]

af = np.linspace(-af_len/2, af_len/2, af_len)*delta_af
mf = np.linspace(-mf_len/2, mf_len/2, mf_len)*delta_mf

plt.subplot(2,2,3)
plt.pcolormesh(af, mf, c, cmap='Blues')
#plt.colorbar()
plt.xlim([af[0], af[-1]])
plt.ylim([mf[0], mf[-1]])
plt.axvline(x=0, color='r', ls='--')
plt.axhline(y=0, color='r', ls='--')
plt.xlabel('$\Delta (a_f/M_f)$')
plt.ylabel('$\Delta (M_f)$')


plt.tight_layout()

plt.show()
