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
      (t_geom, h22r_geom, h22i_geom) = np.loadtxt(data_file, skiprows=1, unpack=True)
    except IOError:
      (t_geom, h22r_geom, h22i_geom) = np.loadtxt(data_file+'.gz', skiprows=1, unpack=True)
    sys.stderr.write('%.3fs.\n'%(time.time() - t_start))
    self.t_geom = t_geom
    self.h22r_geom = h22r_geom
    self.h22i_geom = h22i_geom
    self.data_file = data_file
  def TDwaveform(self, M=20., dist=1000., iota=0., f_low=20., srate=4096., phi_ref=0., f_ref=None, taper=False, trim_right=False, trim_tolerance=1e-6):
    if f_ref is None:
      f_ref = f_low
    t_geom = self.t_geom
    h22r_geom = self.h22r_geom
    h22i_geom = self.h22i_geom
    # Convert from geometric to SI units
    M_geom = M * MTSUN_SI
    t_SI = t_geom * M_geom
    h22r_SI = h22r_geom * M_geom / (dist * 1e6 * PC_SI / C_SI)
    h22i_SI = h22i_geom * M_geom / (dist * 1e6 * PC_SI / C_SI)
    # Compute the complex h22
    h22_SI = h22r_SI - 1j*h22i_SI
    # Compute the phase and take derivative to get the frequency
    phi_SI = np.unwrap(np.angle(h22_SI))
    Foft_SI = np.gradient(phi_SI)/np.gradient(t_SI)/(2*np.pi)
    # Calculate phase at f_ref
    phi_from_F = interpolate.interp1d(Foft_SI, phi_SI)
    phi_at_f_ref = phi_from_F(f_ref)
    phi_at_f_start = phi_SI[0]
    # Correct for phase
    phi0 = np.mod(phi_ref + phi_at_f_start - phi_at_f_ref, 2.*np.pi)
    # compute Y2p2(iota, phi0) and Y2m2(iota, phi0)
    Y2p2 = spinm2_sphharm(2, 2, iota, phi0)
    Y2m2 = spinm2_sphharm(2, -2, iota, phi0)
    # compute the (complex) h(iota, phi0) = h+ - i hx
    hh_SI = h22_SI*Y2p2 + np.conj(h22_SI)*Y2m2
    # get hp and hc
    hp_SI = np.real(hh_SI)
    hc_SI = np.imag(hh_SI)
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

  m1 = 102.797800
  m2 = 136.583600
  M = m1 + m2 # M_SUN
  dist = 5280. # Mpc
  iota = 0.
  phi_ref = 0. # phase at f_ref
  f_ref = 10. # Hz 
  f_low = 10. # Hz
  srate = 4096. # Hz

  modGR_q1_fac1 = modGR_waveform('/home/abhirup/Documents/Work/testGR_IR/src/IHES_EOB_mod/waveforms/popsynth_GR/injection_2739_q_1.33_GR.dat')
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
