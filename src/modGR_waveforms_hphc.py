#!/usr/bin/env python

import time, os, sys, numpy as np
from scipy import interpolate
from lal import MSUN_SI, MTSUN_SI, PC_SI, C_SI
import lal, lalsimulation as lalsim
from spinm2SphHarm import spinm2_sphharm

class modGR_waveform(object):
  def __init__(self, data_file):
    # Load the data
    sys.stderr.write('Loading data .. ')
    t_start = time.time()
    try:
      (t, hp, hc) = np.loadtxt(data_file, skiprows=1, unpack=True)
    except:
      print 'data file not present'
    self.t = t
    self.hp = hp
    self.hc = hc
    self.data_file = data_file
  def TDwaveform(self, M=20., dist=1000., iota=0., f_low=20., srate=4096., phi_ref=0., f_ref=None, taper=False, trim_right=False, trim_tolerance=1e-6):
    if f_ref is None:
      f_ref = f_low
    t_SI = self.t
    hp_SI = self.hp
    hc_SI = self.hc
    # Compute the phase and take derivative to get the frequency
    phi_SI = np.unwrap(np.arctan2(hc_SI, hp_SI))
    Foft_SI = np.gradient(phi_SI)/np.gradient(t_SI)/(2*np.pi)
    hh_SI = hp_SI - 1j*hc_SI
    # Find the frequencies in band
    idx_SI, = np.where(Foft_SI>f_low)
    # trim if asked to
    idx_SI_max = -1
    if trim_right:
      idx_SI_max = np.where(np.abs(hh_SI)/np.max(np.abs(hh_SI))>trim_tolerance)[0][-1]
      #print idx_SI_max
    # create an array of times
    t_arr = np.arange(t_SI[idx_SI[0]], t_SI[idx_SI_max], 1./srate)
    # create an interpolation
    hp_interp = interpolate.interp1d(t_SI, hp_SI)
    hc_interp = interpolate.interp1d(t_SI, hc_SI)
    # arrays to return
    hp_arr = hp_interp(t_arr)
    hc_arr = hc_interp(t_arr)
    # create the REAL8TimeSeries
    hp = lal.CreateREAL8TimeSeries('hp_modGR', 0, f_low, 1./srate, '', len(hp_arr))
    hc = lal.CreateREAL8TimeSeries('hc_modGR', 0, f_low, 1./srate, '', len(hp_arr))
    hp.data.data = hp_arr
    hc.data.data = hc_arr
    # taper waveforms Eq. (3.35) of gr-qc/0001023
    if taper:
      lalsim.SimInspiralREAL8WaveTaper(hp.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
      lalsim.SimInspiralREAL8WaveTaper(hc.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
    return hp, hc
  
if __name__ == '__main__':
  
  import waveforms as wf

  m1 = 36
  m2 = 29
  M = m1 + m2 # M_SUN
  dist = 410. # Mpc
  iota = 0.
  phi_ref = 0. # phase at f_ref
  f_ref = 10. # Hz 
  f_low = 10. # Hz
  srate = 8192. # Hz

  modGR_q1_fac1 = modGR_waveform('../scripts/EccentricTDIMRv1_m1_36_m2_29_f_min_10_emin_0_d_410.dat')
  (t_arr, hp_arr, hc_arr) = wf.data_from_TD(modGR_q1_fac1.TDwaveform(M, dist=dist, iota=iota, f_low=f_low, srate=srate, phi_ref=phi_ref, f_ref=f_ref, taper=True))
  ampl_h_arr = np.sqrt(np.abs(hp_arr)**2.+np.abs(hc_arr)**2.)
  phase_h_arr = np.unwrap(np.angle(hp_arr+1j*hc_arr))

  (t_arr_lalsim, hp_arr_lalsim, hc_arr_lalsim) = wf.data_from_TD(wf.TDwaveform(m1, m2, dist=dist, incl=iota, phi_ref=phi_ref, f_ref=f_ref, f_low=f_low, srate=srate, taper=True))
  ampl_h_arr_lalsim = np.sqrt(np.abs(hp_arr_lalsim)**2.+np.abs(hc_arr_lalsim)**2.)
  phase_h_arr_lalsim = np.unwrap(np.angle(hp_arr_lalsim+1j*hc_arr_lalsim))

  # Comparison of the two
  import matplotlib.pyplot as plt

  plt.figure()
  plt.plot(t_arr-t_arr[0], hp_arr, 'b-', label=r'$h_+$ IHES')
  plt.plot(t_arr-t_arr[0], hc_arr, 'r-', label=r'$h_{\times}$ IHES')
  plt.plot(t_arr-t_arr[0], ampl_h_arr, 'g-', label=r'|h| IHES')
  plt.xlabel('Time (s)')
  plt.legend(loc='upper left')

  plt.figure()
  plt.plot(t_arr_lalsim-t_arr_lalsim[0], hp_arr_lalsim, 'b-', label=r'$h_+$ LALSim')
  plt.plot(t_arr_lalsim-t_arr_lalsim[0], hc_arr_lalsim, 'r-', label=r'$h_{\times}$ LALSim')
  plt.plot(t_arr_lalsim-t_arr_lalsim[0], ampl_h_arr_lalsim, 'g-', label=r'|h| LALSim')
  plt.xlabel('Time (s)')
  plt.legend(loc='upper left')

  plt.figure()
  plt.plot(t_arr-t_arr[0], hp_arr, 'b-', label=r'$h_+$ IHES')
  plt.plot(t_arr_lalsim-t_arr_lalsim[0], hp_arr_lalsim, 'k--', label=r'$h_+$ LALSim')
  plt.xlabel('Time (s)')
  plt.legend(loc='upper left')

  plt.figure()
  plt.plot(t_arr-t_arr[0], hc_arr, 'r-', label=r'$h_{\times}$ IHES')
  plt.plot(t_arr_lalsim-t_arr_lalsim[0], hc_arr_lalsim, 'k--', label=r'$h_\times$ LALSim')
  plt.xlabel('Time (s)')
  plt.legend(loc='upper left')

  plt.figure()
  plt.plot(t_arr-t_arr[0], ampl_h_arr, 'g-', label=r'$|h|$ IHES')
  plt.plot(t_arr_lalsim-t_arr_lalsim[0], ampl_h_arr_lalsim, 'k--', label=r'$|h|$ LALSim')
  plt.xlabel('Time (s)')
  plt.legend(loc='upper left')

  plt.figure()
  plt.plot(t_arr-t_arr[0], phase_h_arr, 'm-', label=r'phase($h$) IHES')
  plt.plot(t_arr_lalsim-t_arr_lalsim[0], phase_h_arr_lalsim, 'k--', label=r'phase($h$) LALSim')
  plt.xlabel('Time (s)')
  plt.legend(loc='upper left')

  min_len = min(len(t_arr), len(t_arr_lalsim))

  plt.figure()
  plt.plot(t_arr_lalsim[:min_len], phase_h_arr[:min_len] - phase_h_arr_lalsim[:min_len], 'm--', label=r'phase($h$) IHES $-$ LALSim')
  plt.xlabel('Time (s)')
  plt.legend(loc='lower left')


  print 'IHES ampl / LALSim ampl = %.2f'%(np.mean(ampl_h_arr[:min_len] / ampl_h_arr_lalsim[:min_len]))

  print 'IHES phase / LALSim phase = %.2f'%(np.mean(phase_h_arr[:min_len] / phase_h_arr_lalsim[:min_len]))

  print 'IHES phase - LALSim phase = %.2f'%(np.mean(phase_h_arr[:min_len] - phase_h_arr_lalsim[:min_len]))

  plt.show()
