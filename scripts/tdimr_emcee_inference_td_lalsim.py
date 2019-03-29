import os
from numpy import sqrt, sin, cos, pi
import matplotlib
matplotlib.use("Pdf")
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import emcee
import sys
sys.path.insert(0, '/home/abhirup/Documents/Work/testGR_IR/src')
import template_lalsiminspiral as phhsi_lal
from pycbc  import  detector
from lal import MSUN_SI, MTSUN_SI, PC_SI, PI, PC_SI, C_SI, GAMMA, MRSUN_SI
import corner
from optparse import OptionParser
import time
import pycbc
import scipy
from scipy import interpolate

def whiten_losc(strain, interp_psd, dt, t0):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)*np.exp(-2*pi*1j*freqs*t0)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_hf[np.isnan(white_hf)] = 0.
    white_hf[np.isinf(white_hf)] = 0.
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

def lnlike(param_vec, data_white, t_data, f_low, dt):
        """
        compute the log likelihood
        
        inputs: 
        param_vec : vector of parameters 
        dr, di, 
        freq : Fourier freq 
        psd : psd vector 
        flow,fcut
        
        output: 
        log_likelhood 
        """
        
	# unpacking the parameter vector 
        Mc, q, dL, i, t0, phi0,  ra, sin_dec, pol= param_vec

        # generate the waveform 
	t_signal, hpt, hct = phhsi_lal.lalsiminspiraltd_waveform_SI("IMRPhenomPv2",Mc,q,dL,i,(phi0 %(2.*pi)),f_low,dt)
	t_signal = np.around(t_signal*1e4)/1e4

        # compute antenna patterns 
        Fp,Fc = detector.overhead_antenna_pattern(ra, np.arcsin(sin_dec), pol)

        signal = Fp*hpt+Fc*hct
	signal_white = whiten_losc(signal, interp_psd, dt, t0)

	t_start = max(t_data[0], t_signal[0])
	t_end = min(t_data[-1], t_signal[-1])

	if cbcstage == 'imr':
	  band_idx_data, = np.where((t_data > t_start) & (t_data < t_end))
	  band_idx_signal, = np.where((t_signal > t_start) & (t_signal < t_end))

	elif cbcstage == 'insp':
	  band_idx_data, = np.where((t_data > t_start) & (t_data <= t_cut))
          band_idx_signal, = np.where((t_signal > t_start) & (t_signal <= t_cut))
	
	elif cbcstage == 'ring':
	  band_idx_data, = np.where((t_data >= t_cut) & (t_data < t_end))
          band_idx_signal, = np.where((t_signal >= t_cut) & (t_signal < t_end))

        like = -np.dot(data_white[band_idx_data]-signal_white[band_idx_signal],data_white[band_idx_data]-signal_white[band_idx_signal])
	return like#log-likelihood


def lnprior(param_vec):
        Mc, q, dL, i, t0, phi0, ra, sin_dec, pol = param_vec
        if 1. < Mc < 200. and 0.05 < q <= 1. and 1.<dL<10000 and 0.<= i <= pi and -15. <= t0 <= 15. and -pi <= phi0 <= 3.*pi and 0. <= ra < 2.*pi and -1. <= sin_dec <= 1. and 0. <= pol <= pi:
                return 2.*np.log(dL)+np.log(np.sin(i))
        return -np.inf



def lnprob(param_vec):
        lp = lnprior(param_vec)
        if not np.isfinite(lp):
                return -np.inf
        return lp + lnlike(param_vec, data_white, t_data, f_low, dt)


##########################################################
###################### MAIN ##############################
##########################################################

start_time = time.time()

# -------------------- inputs -------------------------- # 
parser = OptionParser()
parser.add_option("-d", "--data-fname", dest="data_fname", help="data filename")
parser.add_option("-o", "--out-dir", dest="out_dir", help="output directory")
#parser.add_option("-i", "--init-cond", dest="result", help="initial conditions")
parser.add_option("--t_cut", dest="t_cut", type='float', default=None, help="cutoff time between inspiral and post-inspiral")
parser.add_option("--cbcstage", dest="cbcstage", help="[imr, insp, ring]")
(options, args) = parser.parse_args()
data_fname = options.data_fname
out_dir = options.out_dir
#result = options.result
t_cut = options.t_cut
cbcstage = options.cbcstage
#result = map(float, result.strip('[]').split(','))

os.system('mkdir -p %s'%out_dir)
os.system('cp -r %s %s'%(data_fname, out_dir))
os.system('cp %s %s' %(__file__, out_dir))

f_low = 20.
ndim, nwalkers = 9, 100
num_threads = 30
num_iter = 20000
# ------------------------------------------------------ # 


# read the detector data in Fourier domain. [fourier freq, real part of the data, imaginary part of the data, psd]
t_data, data, hpt, hct = np.loadtxt(data_fname, unpack=True)
t_data = np.array(t_data)
t_data = np.around(t_data*1e4)/1e4
data = np.array(data)
dt = t_data[1] - t_data[0]
f, asd = np.loadtxt('/home/abhirup/Documents/Work/testGR_IR/PSD/AdLIGO_ZDHP.txt', unpack=True)
interp_psd = scipy.interpolate.interp1d(f, asd**2., fill_value=0., bounds_error=False)
data_white = whiten_losc(data, interp_psd, dt, 0.)

debug = True

print '... read and whitened data'

# create initial walkers

result = 28.09555579546043,0.8055555555555556, 500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
mc_init, q_init, dL_init, iota_init, t0_init, phi0_init, ra_init, sin_dec_init, pol_init = result

plt.figure()
plt.plot(t_data, data_white)
plt.savefig(out_dir + '/data_white.png')
plt.close()

pos = [result + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

print '... generated initial walkers. starting sampling...'

# sample the likelihood using EMCEE 
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=num_threads)
sampler.run_mcmc(pos, num_iter)

mc_chain, q_chain, dL_chain, iota_chain, t0_chain, phi0_chain, ra_chain, sin_dec_chain, pol_chain = sampler.chain[:, :, 0].T, sampler.chain[:, :, 1].T, sampler.chain[:, :, 2].T, sampler.chain[:, :, 3].T, sampler.chain[:, :, 4].T, sampler.chain[:, :, 5].T, sampler.chain[:, :, 6].T, sampler.chain[:, :, 7].T, sampler.chain[:, :, 8].T

samples = sampler.chain[:, :, :].reshape((-1, ndim))

#################################################################
# plotting and saving data
#################################################################

# save the data
np.savetxt(out_dir+'/emcee_samples.dat', samples, header='mc q dL i t0 phi0 ra sin(dec) pol')
np.savetxt(out_dir+'/emcee_samples_lnprob.dat', sampler.lnprobability, header='lnprob')

# Inspiral Chain plot
plt.figure(figsize=(16,8))
plt.subplot(621)
plt.plot(mc_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(mc_init + np.std(mc_chain, axis=1), 'r')
plt.axhline(y=mc_init, color='g')
plt.ylabel('mc')
plt.subplot(622)
plt.plot(q_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(q_init + np.std(q_chain, axis=1), 'r')
plt.axhline(y=q_init, color='g')
plt.ylabel('q')
plt.subplot(623)
plt.plot(dL_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(dL_init + np.std(dL_chain, axis=1), 'r')
plt.axhline(y=dL_init, color='g')
plt.ylabel('dL')
plt.subplot(624)
plt.plot(iota_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(iota_init + np.std(iota_chain, axis=1), 'r')
plt.axhline(y=iota_init, color='g')
plt.ylabel('iota')
plt.subplot(625)
plt.plot(t0_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(t0_init + np.std(t0_chain, axis=1), 'r')
plt.axhline(y=t0_init, color='g')
plt.ylabel('t0')
plt.subplot(626)
plt.plot(phi0_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(phi0_init + np.std(phi0_chain, axis=1), 'r')
plt.axhline(y=phi0_init, color='g')
plt.ylabel('phi0')
plt.subplot(627)
plt.plot(ra_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(ra_init + np.std(ra_chain, axis=1), 'r')
plt.axhline(y=ra_init, color='g')
plt.ylabel('ra')
plt.subplot(6,2,8)
plt.plot(sin_dec_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(sin_dec_init + np.std(sin_dec_chain, axis=1), 'r')
plt.axhline(y=sin_dec_init, color='g')
plt.ylabel('dec')
plt.subplot(6,2,9)
plt.plot(pol_chain, color="k", alpha=0.4, lw=0.5)
plt.plot(pol_init + np.std(pol_chain, axis=1), 'r')
plt.axhline(y=pol_init, color='g')
plt.ylabel('pol')
plt.savefig(out_dir + '/samples_chain.png', dpi=300)

# corner plots
plt.figure()
corner.corner(samples, labels=['mc', 'q', 'dL', 'i', 't0', 'phi0', 'ra', 'sin(dec)', 'pol'])
plt.savefig("%s/corner_plot_wo_burnin.png"%out_dir)
plt.close()

# likelihood plots
if debug == True:
    plt.figure(figsize=(16,8))
    plt.subplot(331)
    plt.plot(np.linspace(1,200, 100), [lnlike([mc, q_init, dL_init, iota_init, t0_init, phi0_init, ra_init, sin_dec_init, pol_init], data_white, t_data, f_low, dt) for mc in np.linspace(1, 200, 100)])
    plt.axvline(x=mc_init, color='r')
    plt.subplot(332)
    plt.plot(np.linspace(0.05,1, 100), [lnlike([mc_init, q, dL_init, iota_init, t0_init, phi0_init, ra_init, sin_dec_init, pol_init], data_white, t_data, f_low, dt) for q in np.linspace(0.05, 1, 100)])
    plt.axvline(x=q_init, color='r')
    plt.subplot(333)
    plt.plot(np.linspace(1,10000, 100), [lnlike([mc_init, q_init, dL, iota_init, t0_init, phi0_init, ra_init, sin_dec_init, pol_init], data_white, t_data, f_low, dt) for dL in np.linspace(1, 10000, 100)])
    plt.axvline(x=dL_init, color='r')
    plt.subplot(334)
    plt.plot(np.linspace(0,pi, 100), [lnlike([mc_init, q_init, dL_init, iota, t0_init, phi0_init, ra_init, sin_dec_init, pol_init], data_white, t_data, f_low, dt) for iota in np.linspace(0, pi, 100)])
    plt.axvline(x=iota_init, color='r')
    plt.subplot(335)
    plt.plot(np.linspace(-15,15, 100), [lnlike([mc_init, q_init, dL_init, iota_init, t0, phi0_init, ra_init, sin_dec_init, pol_init], data_white, t_data, f_low, dt) for t0 in np.linspace(-15, 15, 100)])
    plt.axvline(x=t0_init, color='r')
    plt.subplot(336)
    plt.plot(np.linspace(-pi,3.*pi, 100), [lnlike([mc_init, q_init, dL_init, iota_init, t0_init, phi0, ra_init, sin_dec_init, pol_init], data_white, t_data, f_low, dt) for phi0 in np.linspace(-pi, 3*pi, 100)])
    plt.axvline(x=phi0_init, color='r')
    plt.subplot(337)
    plt.plot(np.linspace(0,2*pi, 100), [lnlike([mc_init, q_init, dL_init, iota_init, t0_init, phi0_init, ra, sin_dec_init, pol_init], data_white, t_data, f_low, dt) for ra in np.linspace(0, 2*pi, 100)])
    plt.axvline(x=ra_init, color='r')
    plt.subplot(338)
    plt.plot(np.linspace(-1,1, 100), [lnlike([mc_init, q_init, dL_init, iota_init, t0_init, phi0_init, ra_init, sin_dec, pol_init], data_white, t_data, f_low, dt) for sin_dec in np.linspace(-1, 1, 100)])
    plt.axvline(x=sin_dec_init, color='r')
    plt.subplot(338)
    plt.plot(np.linspace(-1,1, 100), [lnlike([mc_init, q_init, dL_init, iota_init, t0_init, phi0_init, ra_init, sin_dec, pol_init], data_white, t_data, f_low, dt) for sin_dec in np.linspace(-1, 1, 100)])
    plt.axvline(x=sin_dec_init, color='r')
    plt.subplot(339)
    plt.plot(np.linspace(0,pi, 100), [lnlike([mc_init, q_init, dL_init, iota_init, t0_init, phi0_init, ra_init, sin_dec_init, pol], data_white, t_data, f_low, dt) for pol in np.linspace(0, pi, 100)])
    plt.axvline(x=pol_init, color='r')
    plt.savefig(out_dir + '/likelihood.png')


print '... plotted corner plot'

end_time = time.time()
print '... time taken: %.2f seconds'%(end_time-start_time)
