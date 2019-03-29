"""
To generate Adv-LIGO PSD and noise strain.

(C) Walter Del Pozzo. Modified by Archisman Ghosh, 2014-09-18.
"""

import os, glob, numpy as np
import lalinspiral.sbank.psds as psds
from scipy import random, interpolate

lal_prefix = os.getenv('LAL_PREFIX')

def mydiff(f):
  return np.append(np.append([f[1]-f[0]],(f[2:]-f[:-2])/2.),[f[-1]-f[-2]])
  #return 1.

# ### FOURIER TRANSFORMS ###

def fd_from_td(t_arr, strain_t):
  srate = 1./np.mean(np.diff(t_arr))
  # zero-pad to the required length 
  N = int(2**np.ceil(np.log2(len(t_arr))))
  strain_t = np.resize(strain_t, N)
  strain_t[len(t_arr):] = 0
  strain_f = np.fft.rfft(strain_t)/srate
  ff = np.fft.rfftfreq(len(strain_t), t_arr[1]-t_arr[0])
  return (ff, strain_f)

def td_from_fd(f_arr, strain_f):
  N = 2*(len(f_arr)-1)
  srate = N*np.mean(np.diff(f_arr))
  tt = -np.arange(0, N)[::-1]/srate
  strain_t = np.fft.irfft(strain_f)*srate
  return(tt, strain_t)

# ### PSD ###

def LIGOPSD(f, noise_model='aLIGOZeroDetHighPower', nbins=None, is_asd_file=False):
  """
  Returns PSD for an array of frequencies given noise noise_model
  
  Parameters
  ----------
  f : array of frequencies (Hz)
  noise_model : LIGO / Virgo noise model, default = 'aLIGOZeroDetHighPower'
  
  noise_model can be provided in one of the three following ways:
  1. a string recognized by the dictionary lalinspiral.sbank.psds.noise_models
  2. a prefix of noise models in lalsimulation/src like 'aLIGO_DESIGN'
  3. a text file with (f, psd) pairs to interpolate. if (f, asd) pairs present, set is_asd_file=True.
  
  nbins : nbins to smoothen data in PSD file
  is_asd_file : set to True if (f, asd) pairs in file
  
  Returns
  -------
  Array of PSDs
  """
  if np.shape(np.shape(noise_model)) == (2,):
    if nbins is None:
      nbins = len(f)
    (f_data_raw, psd_data_raw) = noise_model
    if is_asd_file:
      psd_data_raw = psd_data_raw**2
    f_data = np.average(f_data_raw[len(f_data_raw)%nbins:].reshape(nbins,-1), axis=1)
    psd_data = np.average(psd_data_raw[len(psd_data_raw)%nbins:].reshape(nbins,-1), axis=1)
    psd_interp = interpolate.interp1d(f_data, psd_data, bounds_error=False)
    return psd_interp(f)
  try:
    return psds.noise_models[noise_model](f)
  except KeyError:
    if os.path.isfile(noise_model):
      noise_file = noise_model
      if is_asd_file:
        (f_data, sqrt_psd_data) = np.loadtxt(noise_file, unpack=True)
        psd_data = sqrt_psd_data**2.
      else:
        (f_data, psd_data) = np.loadtxt(noise_file, unpack=True)
      psd_interp = interpolate.interp1d(f_data, psd_data, bounds_error=False)
    else:
      noise_file = glob.glob(os.path.join(lal_prefix, 'share/lalsimulation/*-%s.txt'%(noise_model)))[0]
      (f_data, sqrt_psd_data) = np.loadtxt(noise_file, unpack=True)
      psd_data = sqrt_psd_data**2.
      psd_interp = interpolate.interp1d(f_data, psd_data, bounds_error=False)
    return psd_interp(f)

# ### Noise realization in FD ###

def FDnoise(f, noise_model='aLIGOZeroDetHighPower', f_low=10., seed=None):
  """
  Returns a noise realization given an array of frequencies given noise noise_model
  
  Parameters
  ----------
  f : array of frequencies (Hz)
  noise_model : LIGO / Virgo noise model, default = 'aLIGOZeroDetHighPower'
  f_low : low frequency cut-off (Hz), default = 0 Hz
  seed : seed for scipy.random random number generator
  
  Returns
  -------
  A complex array containing a FD noise realization
  """
  df = mydiff(f)
  if seed is not None:
    random.seed(seed)
  noise_ampl = 0.5*np.sqrt(LIGOPSD(f, noise_model)/df)
  noise_real = np.array([random.normal(0, noise_ampl[ii]) for ii in range(len(f))])
  noise_imag = np.array([random.normal(0, noise_ampl[ii]) for ii in range(len(f))])
  noise = np.vectorize(complex)(noise_real, noise_imag)
  noise[f<=f_low] = 0.
  return noise

# ### Noise realization in TD ###

def TDnoise(n, dt=1., noise_model='aLIGOZeroDetHighPower', f_low=10., seed=None):
  """
  Returns a time-domain noise realization given properties of the time series given noise noise_model
  
  Parameters
  ----------
  n : length of time series
  dt : real time ealpsed between two adjacent points in the time series (s)
  noise_model : LIGO / Virgo noise model, default = 'aLIGOZeroDetHighPower'
  f_low : low frequency cut-off (Hz), default = 10 Hz
  seed : seed for scipy.random random number generator
  
  Returns
  -------
  Array containing the noise realization in the time domain
  """
  f = np.fft.rfftfreq(int(n), dt)
  nf = FDnoise(f, noise_model=noise_model, f_low=f_low, seed=seed)
  (t, nt) = td_from_fd(f, nf)
  return nt

# ### Noise weighted inner product ###

def inner_product(f, gf, hf, noise_model = 'aLIGOZeroDetHighPower', f_low=10.):
  df = np.mean(np.diff(f))
  band_idx = f>=f_low
  return 4.*(np.sum(np.nan_to_num(df*(np.conj(gf[band_idx])*hf[band_idx])/LIGOPSD(f, noise_model=noise_model)[band_idx])))

# ### SNR calculation for FD waveform ###

def SNR_from_fd(f, hf, noise_model = 'aLIGOZeroDetHighPower', f_low=10., f_high=2048., **kwargs):
  """
  Parameters
  ----------
  f : list of frequencies
  h : FD strain at corresponding frequencies
  noise_model : LIGO / Virgo noise model, default = 'aLIGOZeroDetHighPower'
  f_low : low frequency cut-off (Hz), default = 10 Hz
  
  Returns
  -------
  Optimal SNR
  """
  df = np.mean(np.diff(f))
  band_idx = (f_low<=f) * (f<=f_high)
  return 2*np.sqrt(np.sum(np.nan_to_num(df*abs(hf[band_idx])**2/LIGOPSD(f, noise_model=noise_model, **kwargs)[band_idx])))

# ### PSD for noise ###

def PSD_from_fd(f, hf, nbins=None):
  df = mydiff(f)
  if nbins is None:
    nbins = len(f)
  ff = np.average(f[len(f)%nbins:].reshape(nbins,-1), axis=1)
  psd = 2.*np.average((df*np.abs(hf)**2)[len(hf)%nbins:].reshape(nbins,-1), axis=1)
  return (ff, psd)

""" Usage """
if __name__=="__main__":
  import noiseutils as nu
  from scipy.signal import welch
  import matplotlib.pyplot as plt
  try:
    from pycbc.types import TimeSeries
    from pycbc.psd import estimate
  except ImportError:
    print "Unable to import PyCBC dependencies. Will proceed without PyCBC."
    pycbc_flag = False
  else:
    pycbc_flag = True
  
  # Plot a few noise curves
  ff = np.linspace(8, 4096, 32768/4)
  noises = {'aLIGOZeroDetHighPower': 'ko-', 'aLIGO_EARLY_HIGH': 'r:', 'aLIGO_MID_HIGH': 'r--', 'aLIGO_LATE_HIGH': 'r-', 'aLIGO_DESIGN': 'r.-', 'AdV_EARLY_HIGH': 'g:', 'AdV_MID_HIGH': 'g--', 'AdV_LATE_HIGH': 'g-', 'AdV_DESIGN': 'g.-'}
  plt.figure()
  for (noise_model, fmt_string) in noises.items():
    plt.loglog(ff, LIGOPSD(ff, noise_model), fmt_string, label=noise_model)
  plt.legend(loc='best')
  plt.xlim(1e1, 1e4)
  plt.xlabel(r'Frequency (Hz)')
  plt.ylabel(r'PSD (Hz$^{-1}$)')
  plt.title('Noise curves')
  
  # Plot a noise realization and re-estimate its PSD
  noise_model = 'aLIGOZeroDetHighPower'
  seglen=8
  ff = np.arange(16, 4096, 1./seglen)
  initial_psd = nu.LIGOPSD(ff, noise_model)
  noise_fd = nu.FDnoise(ff, noise_model)
  (fbinned, reestimated_psd) = nu.PSD_from_fd(ff, noise_fd, nbins=1000)
  (tt, noise_td) = nu.td_from_fd(ff, noise_fd)
  dt = tt[1]-tt[0]
  if pycbc_flag:
    w_est_pycbc = estimate.welch(TimeSeries(noise_td, dt))
  (f_scipy, w_est_scipy) = welch(noise_td, 1./dt)
  plt.figure()
  plt.loglog(ff, noise_fd, 'g.', label='noise realization with FDnoise')
  plt.loglog(fbinned, np.sqrt(reestimated_psd), 'r--', linewidth=2, label='re-estimated PSD with PSD_from_fd')
  if pycbc_flag:
    plt.loglog(w_est_pycbc.get_sample_frequencies(), np.sqrt(w_est_pycbc), 'c-', label='re-estimated PSD with pycbc.psd.estimate.welch')
  plt.loglog(f_scipy, np.sqrt(w_est_scipy), 'b-', label='re-estimated PSD with scipy.signal.welch')
  plt.loglog(ff, np.sqrt(initial_psd), 'k-', linewidth=2, label=noise_model)
  plt.legend(loc='best')
  plt.xlim(1e1, 1e4)
  plt.xlabel(r'Frequency (Hz)')
  plt.ylabel(r'ASD (Hz$^{-1/2}$)')
  plt.title('Noise realization and re-estimation of PSD')
    
  # Calculate SNR for a LALSimulation waveform
  import waveforms as wf
  m1 = 50.
  m2 = 50.
  dist = 1000.
  iota = 0.
  ra = 0.
  dec = 0.
  psi = 0.
  f_low = 10.
  f_start = 20.
  srate = 4096.
  f_stop = srate/2.
  trig_time = 1126285216
  approx = 'EOBNRv2'
  plt.figure()
  noise_model = 'aLIGO_EARLY_HIGH'
  plt.loglog(ff, np.sqrt(LIGOPSD(ff, noise_model)), 'k--', linewidth=2, label=noise_model)
  noise_model = 'aLIGOZeroDetHighPower'
  plt.loglog(ff, np.sqrt(LIGOPSD(ff, noise_model)), 'k-', linewidth=2, label=noise_model)
  for det in ['H1', 'L1']:
    (ff, strain_f) = fd_from_td(*wf.detector_strain(det, ra, dec, psi, trig_time, *wf.data_from_TD(wf.TDwaveform(m1, m2, dist=dist, incl=iota, approx=approx, f_low=f_low, srate=srate))))
    noise_model = 'aLIGO_EARLY_HIGH'
    aleh_snr = SNR_from_fd(ff, strain_f, noise_model, f_low=f_start)
    noise_model = 'aLIGOZeroDetHighPower'
    zdhp_snr = SNR_from_fd(ff, strain_f, noise_model, f_low=f_start)
    plt.loglog(ff, 2*np.abs(strain_f)*np.sqrt(ff), '-', label=r'%s: ALEH SNR = %.2f, ZDHP SNR = %.2f'%(det, aleh_snr, zdhp_snr))
  plt.legend(loc='lower left')
  plt.xlim(1e1, 1e4)
  plt.xlabel(r'Frequency, $f$ (Hz)')
  plt.ylabel(r'ASD and $2|h(f)|\sqrt{f}$ (Hz$^{-1/2}$)')
  plt.title('SNR for LALSimulation waveform')
  
  plt.show()
