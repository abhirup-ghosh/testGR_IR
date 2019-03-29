"""
Wrapper for LALSimulation waveforms, projection over antenna patterns using pylal utilities.

Archisman Ghosh, Abhirup Ghosh, P. Ajith; 2013-05-29; last modified: 2016-11-03

$Id:$
"""

import sys, math, numpy as np
import lal, lalsimulation as lalsim
from lal import MSUN_SI, PC_SI
from pylal import antenna, inject
from pylal.date import XLALTimeDelayFromEarthCenter
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

# ### EXTRACT DATA FROM SWIG OBJECTS ###

def data_from_TD(h):
  if np.shape(h) == ():
    # define time
    t = -np.arange(0, h.data.length)[::-1]*h.deltaT
    return (t, h.data.data)
  elif np.shape(h) == (2,):
    t = -np.arange(0, h[0].data.length)[::-1]*h[0].deltaT
    return (t, h[0].data.data, h[1].data.data)

def data_from_FD(hf):
  if np.shape(hf) == ():
    # define frequency
    f = np.arange(0., hf.data.length)*hf.deltaF
    return (f, hf.data.data)
  elif np.shape(hf) == (2,):
    f = np.arange(0., hf[0].data.length)*hf[0].deltaF
    return (f, hf[0].data.data, hf[1].data.data)

# ### FOURIER TRANSFORMS ###

def FD_from_TD(h):
  if np.shape(h) == (2,):
    return FD_from_TD(h[0]), FD_from_TD(h[1])
  srate = 1./h.deltaT
  # zero-pad to the required length 
  N = int(2**np.ceil(np.log2(h.data.length)))
  h = lal.ResizeREAL8TimeSeries(h, 0, N)
  # create vector to hold the FFT and FFT plan, do the FFT
  df = srate/N
  hf = lal.CreateCOMPLEX16FrequencySeries("h(f)", h.epoch, h.f0, df, lal.HertzUnit, int(N/2+1))
  fftplan = lal.CreateForwardREAL8FFTPlan(N, 0)
  lal.REAL8TimeFreqFFT(hf, h, fftplan)
  # return
  return hf

def TD_from_FD(hf):
  if np.shape(hf) == (2,):
    return TD_from_FD(hf[0]), TD_from_FD(hf[1])
  N = 2*(hf.data.length-1)
  srate = N*hf.deltaF
  h = lal.CreateREAL8TimeSeries("h(t)", hf.epoch, 0., 1./srate, lal.DimensionlessUnit,N)
  fftplan = lal.CreateReverseREAL8FFTPlan(N, 0)
  lal.REAL8FreqTimeFFT(h, hf, fftplan)
  return h

# ### LALSIMULATION WAVEFORMS ###

def TDwaveform(m1=10., m2=10., s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., dist=1000., incl=0., approx='EOBNRv2', amplO=-1, phaseO=-1, srate=4096., f_low=40., f_high=None, phi_ref=0., f_ref=100., taper=False, nonGR=None):
  approx_enum = lalsim.GetApproximantFromString(approx)
  if f_high is None:
    f_high = 0.5*srate
  # generate template waveforms in the time domain
  (hp, hc) = lalsim.SimInspiralChooseTDWaveform(
    phi_ref,        # phi_ref (rad)
    1./srate,       # delta_t (sec)
    m1 * MSUN_SI,   # mass 1 (kg)
    m2 * MSUN_SI,   # mass 2 (kg)
    s1x,            # spin1x (S1x/m1^2)  
    s1y,            # spin1y (S1y/m1^2)
    s1z,            # spin1z (S1z/m1^2)
    s2x,            # spin2x (S2x/m2^2)
    s2y,            # spin2y (S2y/m2^2)
    s2z,            # spin2z (S2z/m2^2)
    f_low,          # f_low (Hz)
    f_ref,         # f_ref (Hz)
    dist * 1.0e6 * PC_SI,  # distance (m)
    incl,           # incl_angle angle (rad)
    0.,             # lambda1 (Love number)
    0.,             # lambda2 (Love number)
    None,           # waveform flags
    nonGR,           # non-GR parameters
    amplO,          # amplitude PN order
    phaseO,         # phase PN order
    approx_enum     # approximant
    )
  # taper waveforms Eq. (3.35) of gr-qc/0001023
  if taper:
    lalsim.SimInspiralREAL8WaveTaper(hp.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
    lalsim.SimInspiralREAL8WaveTaper(hc.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
  return hp, hc

def FDwaveform(m1=10., m2=10., s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., dist=1000., incl=0., approx='IMRPhenomB', amplO=-1, phaseO=-1, seglen=32., srate=4096., f_low=40., f_high=None, phi_ref=0., f_ref=100., nonGR=None):
  approx_enum = lalsim.GetApproximantFromString(approx)
  if f_high is None:
    f_high = 0.5*srate
  # generate template waveforms in the time domain
  (hp, hc) = lalsim.SimInspiralChooseFDWaveform(
    phi_ref,        # phi_ref (rad)
    1./seglen,      # delta_f (Hz)
    m1 * MSUN_SI,   # mass 1 (kg)
    m2 * MSUN_SI,   # mass 2 (kg)
    s1x,            # spin1x (S1x/m1^2)  
    s1y,            # spin1y (S1y/m1^2)
    s1z,            # spin1z (S1z/m1^2)
    s2x,            # spin2x (S2x/m2^2)
    s2y,            # spin2y (S2y/m2^2)
    s2z,            # spin2z (S2z/m2^2)
    f_low,          # f_low (Hz)
    f_high,         # f_high (Hz)
    f_ref,         # f_ref (Hz)
    dist * 1.0e6 * PC_SI,  # distance (m)
    incl,           # incl_angle angle (rad)
    0.,             # lambda1 (Love number)
    0.,             # lambda2 (Love number)
    None,           # waveform flags
    nonGR,           # non-GR parameters
    amplO,          # amplitude PN order
    phaseO,         # phase PN order
    approx_enum     # approximant
    )
  return hp, hc

# ### PROJECTION OVER ANTENNA PATTERN FUNCTIONS ###

def detector_strain(det, ra, dec, psi, trig_time, t_arr, hp_arr, hc_arr, seglen=None, gps_end=None):
  srate = 1./(t_arr[1] - t_arr[0])
  wf_len = t_arr[-1] - t_arr[0]
  if gps_end is None:
    gps_end = trig_time + 2
  if seglen is None:
    seglen = int(2**math.ceil(np.log2(wf_len + 2)))
  else:
    seglen = int(seglen)
  gps_start = gps_end - seglen
  init_padding = seglen - (gps_end - trig_time) - wf_len
  init_padding_idx = int(init_padding * srate)
  # Calculate the time shift
  timedelay = XLALTimeDelayFromEarthCenter(inject.cached_detector_by_prefix[det].location, ra, dec, LIGOTimeGPS(trig_time))
  timedelay_idx = int(round(timedelay * srate))
  full_padding_idx = init_padding_idx + timedelay_idx
  # Calculate the antenna pattern functions
  Fplus, Fcross, Faverage, Qvalue = antenna.response(trig_time, ra, dec, 0., psi, 'radians', det)
  strain_t = np.zeros(int(seglen*srate))
  for k in xrange(len(strain_t)):
    if full_padding_idx < k < full_padding_idx + len(t_arr):
      strain_t[k] = Fplus*hp_arr[k-full_padding_idx] + Fcross*hc_arr[k-full_padding_idx]
  return (np.arange(0, len(strain_t))/srate, strain_t)

""" Usage """
if __name__=="__main__":
  import waveforms as wf
  import noiseutils as nu
  import matplotlib.pyplot as plt
  
  # Comparison of LAL and numpy FFTs
  
  # FD waveforms
  (f, hf) = wf.data_from_FD(wf.FDwaveform(20., 20., f_low=10.)[0])
  plt.loglog(f, np.abs(hf), label='FD flow 10')
  (f, hf) = wf.data_from_FD(wf.FDwaveform(20., 20., f_low=20.)[0])
  plt.loglog(f, np.abs(hf), label='FD flow 20')
  (f, hf) = wf.data_from_FD(wf.FDwaveform(20., 20., f_low=40.)[0])
  plt.loglog(f, np.abs(hf), label='FD flow 40')

  # LAL FFT of TD waveforms
  (f, hf) = wf.data_from_FD(wf.FD_from_TD(wf.TDwaveform(20., 20., f_low=10.)[0]))
  plt.loglog(f, np.abs(hf), label='LAL FFT flow 10')
  (f, hf) = wf.data_from_FD(wf.FD_from_TD(wf.TDwaveform(20., 20., f_low=20.)[0]))
  plt.loglog(f, np.abs(hf), label='LAL FFT flow 20')
  (f, hf) = wf.data_from_FD(wf.FD_from_TD(wf.TDwaveform(20., 20., f_low=40.)[0]))
  plt.loglog(f, np.abs(hf), label='LAL FFT flow 40')

  # numpy fft of TD waveforms
  (f, hf) = nu.fd_from_td(*wf.data_from_TD(wf.TDwaveform(20., 20., f_low=10.)[0]))
  plt.loglog(f, abs(hf), label='np.fft flow 10')
  (f, hf) = nu.fd_from_td(*wf.data_from_TD(wf.TDwaveform(20., 20., f_low=20.)[0]))
  plt.loglog(f, abs(hf), label='np.fft flow 20')
  (f, hf) = nu.fd_from_td(*wf.data_from_TD(wf.TDwaveform(20., 20., f_low=40.)[0]))
  plt.loglog(f, abs(hf), label='np.fft flow 40')
  
  plt.legend(loc='best')
  plt.xlim(1e1, 1e4)
  plt.xlabel(r'Frequency (Hz)')
  plt.ylabel(r'$h(f)$ (Hz$^{-1}$)')
  plt.title('Comparison of LAL and numpy FFTs')
  
  plt.show()
