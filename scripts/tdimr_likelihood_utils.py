import matplotlib as mpl
mpl.use('Agg')
import commands, os, lal, glob, numpy as np
import matplotlib.pyplot as plt
import math
import standard_gwtransf as gw
from pylal import antenna
import scipy 
from scipy import interpolate, random
import lal
import pycbc.psd
import pycbc.types
import pycbc.filter
import pycbc.noise

G = 6.67408e-11
M_sun = 1.99e30
AU = 1.496e+11
Pc = 3.086e+16
c = 3e8

######################################################################################
# function defining a chirp signal
######################################################################################
def chirp(t, Mc, dL):
	tc = 0.
	phic = 0.
	dL = dL * 1e6 * Pc
	Mc = Mc * M_sun

	tau = -t
	phi = -2. * ((5. * G * Mc / c**3. / tau)**(-5./8.)) + phic
	A = (1./(2.*dL)) * (G * Mc / c**2.)**(5./4.) * (5./(c*tau))**(1./4.) # inclination is pi/2
	h = A * np.cos(phi)
	return h

######################################################################################
# function defining the QNM signal corresponding to the
# dominant (2,2) mode
######################################################################################
def qnm(t, t0_ring, A, tau, f_qnm, phi0):
        return A*np.exp(-(t - t0_ring)/tau)*np.cos(2.*np.pi*f_qnm*(t - t0_ring) - phi0)

######################################################################################
# LALSimulation PSD 
######################################################################################

def psd_from_txt(asd_file='ZERO_DET_high_P', delta_f=1.0/16., flow=10., f_high=1024.):
      
    asd_file = glob.glob('/home/abhirup/Documents/Work/testGR_IR/PSD/*%s*.txt'%(asd_file))[0]
    flen = int(f_high / delta_f) + 1
    psd = pycbc.psd.from_txt(asd_file, flen, delta_f, flow, is_asd_file=True)
    return psd

######################################################################################
# ### Noise realization in TD ###
######################################################################################

def TDnoise(delta_t=1./2048., noise_length=16., asd_file='ZERO_DET_high_P', flow=10., f_high=1024.):
    delta_f = 1./noise_length
    tsamples = int(noise_length / delta_t)
    psd = psd_from_txt(asd_file, delta_f, flow, f_high)
    nt = pycbc.noise.noise_from_psd(tsamples, delta_t, psd, seed=127)
    return nt

######################################################################################
# Estimating the PSD of a time series
######################################################################################

def PSD_from_TDnoise(nt, delta_t):
    # Estimate the PSD
    ## We'll choose 4 seconds PSD samples that are overlapped 50 %
    seg_len = int(4 / delta_t)
    seg_stride = seg_len / 2
    estimated_psd = pycbc.psd.welch(nt, seg_len=seg_len, seg_stride=seg_stride)
    return estimated_psd

######################################################################################
# matched filter SNR
######################################################################################

def SNR_from_td(h, delta_t, psd, flow, srate):
    h = pycbc.types.timeseries.TimeSeries(h, delta_t=delta_t)
    snr = pycbc.filter.matchedfilter.sigma(h, psd=psd, low_frequency_cutoff=flow, high_frequency_cutoff=srate/2.)
    return snr

######################################################################################
# function to whiten data given an interpolated PSD and dt
######################################################################################
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)

    # whitening: transform to freq domain, divide by asd, then transform back, 
    # taking care to get normalization right.
    hf = np.fft.rfft(strain)
    white_hf = hf / (np.sqrt(interp_psd(freqs) /dt/2.))
    white_hf[np.isinf(white_hf)] = 0.
    white_hf[np.isnan(white_hf)] = 0.
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht


######################################################################################
# function defining a basic MCMC module
######################################################################################
def metrohast(m1_old, m2_old, t, d, phase, N_steps):
    m1_samples, m2_samples = [], []
    for idx in range(N_steps):
      t_old, h_old = wavegen_from_lalsiminspiral(approx, m1_old, m2_old, dL, flow, srate)
      h_interp = scipy.interpolate.interp1d(t_old, h_old, fill_value=0., bounds_error=False)
      h_old = h_interp(t)
      l_old = -0.5*np.dot(d - h_old, d - h_old)/(sigma*sigma)
      m1_new, m2_new = m1_old + np.random.rand(), m2_old + np.random.randn()
      if comp_mass_min < m1_new < comp_mass_max and  comp_mass_min < m2_new < comp_mass_max:
        t_new, h_new = wavegen_from_lalsiminspiral(approx, m1_new, m2_new, dL, flow, srate)
        h_interp = scipy.interpolate.interp1d(t_new, h_new, fill_value=0., bounds_error=False)
        h_new = h_interp(t)
        l_new = -0.5*np.dot(d - h_new, d - h_new)/(sigma*sigma)
        R = l_new - l_old
        if np.exp(R) > np.random.random():
                m1_old, m2_old, l_old = m1_new, m2_new, l_new
                print m1_old, m2_old, l_old
                m1_samples, m2_samples = np.append(m1_samples, m1_old), np.append(m2_samples, m2_old)

    return m1_samples, m2_samples

######################################################################################
# function to generate a lalsim-inspiral time domain waveform
# and then the time-domain strain
######################################################################################
def signalTD(approx, m1, m2, dL, phiref, flow, srate, gpsTime, rightAscension, declination, inclination, polarization, unit, detector):
    s = commands.getoutput('lalsim-inspiral --domain="time" --approximant=%s --m1 %f --m2 %f --distance %f --phiRef %f --f-min %f --sample-rate %d --inclination %f'%(approx, m1, m2, dL, phiref, flow, srate, inclination)).split('\n');
    s = s[1:]
    t_h, hp, hc = np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s))
    for idx in range(len(s)):
        t_h[idx], hp[idx], hc[idx] = s[idx].split('\t')
    fp, fc, favg, qvalue = antenna.response(gpsTime, rightAscension, declination, inclination, polarization, unit, detector)
    h = fp*hp + fc*hc
    return t_h, h

#############i#########################################################################
# function to generate a lalsim-inspiral frequency domain waveform
######################################################################################
def signalFD(approx, m1, m2, dL, phiref, flow, srate, gpsTime, rightAscension, declination, inclination, polarization, unit, detector):
    s = commands.getoutput('lalsim-inspiral -F --domain="freq" --approximant=%s --m1 %f --m2 %f --distance %f --phiRef %f --f-min %f --sample-rate %d --inclination %f'%(approx, m1, m2, dL, phiref, flow, srate, inclination)).split('\n');
    s = s[1:]
    f, hp_real, hp_imag, hc_real, hc_imag = np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s))
    for idx in range(len(s)):
        f[idx], hp_real[idx], hp_imag[idx], hc_real[idx], hc_imag[idx] = s[idx].split('\t')
    hp = hp_real + 1j*hp_imag
    hc = hc_real + 1j*hc_imag
    fp, fc, favg, qvalue = antenna.response(gpsTime, rightAscension, declination, inclination, polarization, unit, detector)
    hf = fp*hp + fc*hc
    return f, hf

def wavegen_from_lalsiminspiralFD(approx, m1, m2, a1z, a2z, dL, flow, srate):
    s = commands.getoutput('lalsim-inspiral -F --domain="freq" --approximant=%s --m1 %f --m2 %f --spin1z %f --spin2z %f --distance=%f --f-min %f --sample-rate %d'%(approx, m1, m2, a1z, a2z, dL, flow, srate)).split('\n');
    s = s[1:]
    f, hp_real, hp_imag, hc_real, hc_imag = np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s)), np.zeros(len(s))
    for idx in range(len(s)):
        f[idx], hp_real[idx], hp_imag[idx], hc_real[idx], hc_imag[idx] = s[idx].split('\t')
    hp = hp_real + 1j*hp_imag
    hc = hc_real + 1j*hc_imag
    return f, hp, hc

######################################################################################
# function to generate the ET and adLIGO PSDs
######################################################################################
# Slide 7 of https://dcc.ligo.org/DocDB/0119/G1500703/002/150522-GWADW.pdf (https://arxiv.org/pdf/1410.0612.pdf)

def CE_psd(f_start, srate):
	f_final = srate / 2.
	f = np.arange(f_start, f_final, 1.) 

	Sf = 10.**(-50.)*(11.5*(f/10.)**(-50) + (f/25.)**(-10) + (f/53.)**(-4) + 2.*(f/80.)**(-2.) + 2.*(f/100.)**2.)
	return f, Sf

# Eq. 2.2, 2.3 from https://arxiv.org/pdf/1005.0304.pdf

def ET_psd(f_start, srate):
    if f_start >= 1.:
	
	f_final = srate / 2.	

	a1 = 2.39*1e-27
	a2 = 0.349
	a3 = 1.76
	a4 = 0.409 
	b1 = -15.64
	b2 = -2.145
	b3 = -0.12
	b4 = 1.10
	f0 = 100.
	S0 = 1e-50

	f = np.arange(f_start, f_final, 1.)
	x = f/f0
	Sf = S0 * (a1*(x**b1) + a2*(x**b2) + a3*(x**b3) + a4*(x**b4)) * (a1*(x**b1) + a2*(x**b2) + a3*(x**b3) + a4*(x**b4))
	return f, Sf
    else:
	return np.inf

# Eq. 2.2, 2.3 from https://arxiv.org/pdf/1005.0304.pdf
def adLIGO_psd(f_start, srate):
    if f_start >= 10.:

	f_final = srate/2.
	f = np.arange(f_start, f_final, 1.)

	S0 = 1e-49
	f0 = 215.
	x = f/f0
	Sf = S0 * (10. ** (16. - 4*(f-7.9)**2.) + 2.4 * 1e-62 * x **(-50) + 0.08 * x ** (-4.69) + 123.35 * (1. - 0.23*x**2. + 0.0764 * x**4.) / (1. + 0.17 * x**2.))
	return f, Sf
    else:
	return np.inf

# Eq 6, https://arxiv.org/pdf/1202.4031.pdf
def adVIRGO_asd(f_start, srate):
    if f_start >= 10.:

        f_final = srate/2.
        f = np.arange(f_start, f_final, 1.)

	S0_root = 1.259 * 1e-24
        f0 = 300.
        x = np.log(f/f0)

	Sf_root = S0_root * (0.07*np.exp(-0.142 - 1.437*x + 0.407*x*x) + 3.10*np.exp(-0.466 - 1.043*x - 0.548*x*x) + 0.40*np.exp(-0.304 + 2.896*x - 0.293*x*x) + 0.09*np.exp(1.466 + 3.722*x - 0.984*x*x))
	return f, Sf_root
    else:
	return np.inf

######################################################################################
# function to obtain (Mf, af) from (f_qnm, tau) by inverting qnmfreqs_berti
# (f_qnm, tau) are in SI units (Hz, seconds)
######################################################################################

def mfaf_from_qnmfreqs_berti(f_qnm, tau, l, m, n):

        # load the data file containing the fits (Berti et al,  gr-qc/0512160)
        lVec, mVec, nVec, f1Vec, f2Vec, f3Vec, q1Vec, q2Vec, q3Vec = np.loadtxt('../src/Berti_QNMfitcoeffsWEB.dat', unpack=True)

        idx = np.logical_and(np.logical_and(lVec == l, mVec == m), nVec == n)

        # evaluate the Berti et al fits to the complex frequencies 
        if len(lVec[idx]) == 1:

                 f1 = f1Vec[idx]
                 f2 = f2Vec[idx]
                 f3 = f3Vec[idx]
                 q1 = q1Vec[idx]
                 q2 = q2Vec[idx]
                 q3 = q3Vec[idx]

	omega = 2.*np.pi*f_qnm
	af = 1. - (((omega*tau/2.) - q1)/q2)**(1./q3)
	Mf = (f1 + f2*(1. - af)**f3)/omega

	return Mf/lal.MTSUN_SI, af

######################################################################################
# defining confidence class
######################################################################################
class confidence(object):
  def __init__(self, counts):
    # Sort in descending order in frequency
    self.counts_sorted = np.sort(counts.flatten())[::-1]
    # Get a normalized cumulative distribution from the mode
    self.norm_cumsum_counts_sorted = np.cumsum(self.counts_sorted) / np.sum(counts)
    # Set interpolations between heights, bins and levels
    self._set_interp()
  def _set_interp(self):
    self._length = len(self.counts_sorted)
    # height from index
    self._height_from_idx = interpolate.interp1d(np.arange(self._length), self.counts_sorted, bounds_error=False, fill_value=0.)
    # index from height
    self._idx_from_height = interpolate.interp1d(self.counts_sorted[::-1], np.arange(self._length)[::-1], bounds_error=False, fill_value=self._length)
    # level from index
    self._level_from_idx = interpolate.interp1d(np.arange(self._length), self.norm_cumsum_counts_sorted, bounds_error=False, fill_value=1.)
    # index from level
    self._idx_from_level = interpolate.interp1d(self.norm_cumsum_counts_sorted, np.arange(self._length), bounds_error=False, fill_value=self._length)
  def level_from_height(self, height):
    return self._level_from_idx(self._idx_from_height(height))
  def height_from_level(self, level):
    return self._height_from_idx(self._idx_from_level(level))

######################################################################################
# defining optimal snr module
######################################################################################
def optimal_snr_module(filename):
        data = np.genfromtxt(filename, dtype=None, names=True, usecols=(0,1,2,3,4,5,6))
        var_names = [d[0] for d in data]
        stat_names = data.dtype.names
        optimal_snr = data[var_names.index('optimal_snr')][stat_names.index('mean')+1]
        h1_optimal_snr = data[var_names.index('h1_optimal_snr')][stat_names.index('mean')+1]
        l1_optimal_snr = data[var_names.index('l1_optimal_snr')][stat_names.index('mean')+1]
        v1_optimal_snr = data[var_names.index('v1_optimal_snr')][stat_names.index('mean')+1]
        m1_rec = data[var_names.index('m1')][stat_names.index('mean')+1]
        m2_rec = data[var_names.index('m2')][stat_names.index('mean')+1]
        a1z_rec = data[var_names.index('a1z')][stat_names.index('mean')+1]
        a2z_rec = data[var_names.index('a2z')][stat_names.index('mean')+1]
        Mf_rec = 0.#data[var_names.index('mf')][stat_names.index('mean')+1]
        af_rec = data[var_names.index('af')][stat_names.index('mean')+1]
        return optimal_snr, h1_optimal_snr, l1_optimal_snr, v1_optimal_snr, m1_rec, m2_rec, a1z_rec, a2z_rec, Mf_rec, af_rec
