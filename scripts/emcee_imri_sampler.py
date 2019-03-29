#import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import h5py
import emcee
import corner
import utility_codes as uc
import time
import commands

start_time = time.time()

def wavegen_from_lalsiminspiralFD(approx, m1, m2, dL, flow, fhigh, srate):
    s = commands.getoutput('lalsim-inspiral --frequency-domain --domain="freq" --approximant=%s --m1 %f --m2 %f --distance=%f --f-min %f --sample-rate %d'%(approx, m1, m2, dL, flow, srate)).split('\n');
    s = s[1:]

    f, hpr, hpi, hcr, hci = np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s))

    for idx in range(len(s)):
        f[idx], hpr[idx], hpi[idx], hcr[idx], hci[idx] = s[idx].split('\t')
    hp = hpr - 1j*hpi
    hc = hcr - 1j*hci
    hf = hp - 1j*hc
    return f, hf

def lnlike(theta, f, d):
    m1, m2 = theta
    freqs, hf = wavegen_from_lalsiminspiralFD(approx, m1, m2, dL, flow, srate)
    hf_interp = scipy.interpolate.interp1d(freqs, hf, fill_value=0., bounds_error=False)
    hf = h_interp(f)
    return -0.5*np.dot(d[:-1]-hf[:-1], d[:-1]-hf[:-1] )

def lnprior(theta):
    m1, m2 = theta
    if 1 < m1 < 100 and 1 < m2 < 100:
        return 0.0
    return -np.inf

def lnprob(theta, f, data):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, f, data)



# definition of constant
G = 6.67408e-11
M_sun = 1.99e30
AU = 1.496e+11
Pc = 3.086e+16
c = 3e8

approx = 'IMRPhenomPv2'
srate = 1024
flow = 20.
dL = 1000
m1_inj, m2_inj = 30., 30.

f, hf = wavegen_from_lalsiminspiralFD(approx, m1_inj, m2_inj, dL, flow, srate)
plt.plot(f, hf)
plt.show()
exit()


# whitened data
strain_whiten = uc.whiten(strain,psd,dt)

plt.figure()
plt.plot(time, strain, 'r')
plt.plot(time, strain_whiten * (np.sqrt(1e-48 /dt/2.)), 'g--')

# MCMC
ndim, nwalkers = 3, 500
pos = [[r_inj, m1_inj, m2_inj] + 10.*np.random.rand(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(time[:-1], strain_whiten[:-1]))

sampler.run_mcmc(pos, 5000)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

plt.figure()
corner.corner(samples, labels=["$r$", "$m_1$", "$m_2$"],
                      truths=[r_inj, m1_inj, m2_inj])
plt.savefig('../plots/emcee_imri_postwhitening_walkers500_steps1500.png', dpi=300)

np.savetxt('../plots/samples_emcee_walkers500_steps1500.dat', np.c_[sampler.flatchain[:,0], sampler.flatchain[:,1], sampler.flatchain[:,2]])

end_time = time.time()
print '... completed in %.2f seconds'%(end_time - start_time)
