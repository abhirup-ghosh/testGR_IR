#!/usr/bin/python

"""
Script to generate PN waveforms using the LALSimulation library

P. Ajith, 2013-05-29

$Id:$
"""
import numpy as np
import advLIGOpsd
import lalsimulation as lalsim
from lal import MSUN_SI, PC_SI

def TDwaveform(m1=10., m2=10., s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., dist=1000., incl=0., approx='EOBNRv2', amplO=-1, phaseO=-1, seglen=32., srate=4096., f_low=40., f_high=None, phi_ref=0., taper=False):
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
    f_high,         # f_ref (Hz)
    dist * 1.0e6 * PC_SI,  # distance (m)
    incl,           # incl_angle angle (rad)
    0.,             # lambda1 (Love number)
    0.,             # lambda2 (Love number)
    None,           # waveform flags
    None,           # non-GR parameters
    amplO,          # amplitude PN order
    phaseO,         # phase PN order
    approx_enum     # approximant
    )
  # taper waveforms Eq. (3.35) of gr-qc/0001023
  if taper:
    lalsim.SimInspiralREAL8WaveTaper(hp.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
    lalsim.SimInspiralREAL8WaveTaper(hc.data, lalsim.SIM_INSPIRAL_TAPER_STARTEND)
  # define time
  t = np.linspace(0, hp.data.length/srate, hp.data.length)
  return (t, hp.data.data, hc.data.data)

def FDwaveform(m1=10., m2=10., s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., dist=1000., incl=0., approx='IMRPhenomB', amplO=-1, phaseO=-1, seglen=32., srate=4096., f_low=40., f_high=None, phi_ref=0., f_ref=100.):
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
    None,           # non-GR parameters
    amplO,          # amplitude PN order
    phaseO,         # phase PN order
    approx_enum     # approximant
    )
  # define frequency
  f = np.linspace(0., f_high, hc.data.length)
  return (f, hp.data.data, hc.data.data)

""" Usage """
if __name__=="__main__":
  import matplotlib.pyplot as plt
  
  (tt, hpt, hct) = TDwaveform(10., 10., approx='EOBNRv2', dist=1000.)
  noise_TD = LIGOpsd.TDnoise(len(tt), tt[1]-tt[0], f_low=10.)
  
  fff = np.fft.rfftfreq(len(tt), tt[1]-tt[0])
  hptf = np.fft.rfft(hpt)/(len(tt))
  hctf = np.fft.rfft(hct)/(len(tt))
  
  (ff, hpf, hcf) = FDwaveform(10., 10., approx='IMRPhenomB', dist=1000.)
  
  plt.figure('EOBNRv2 waveform')
  plt.title('EOBNRv2 waveform')
  plt.plot(tt, hpt, 'r-', label='$h_+$')
  plt.plot(tt, hct, 'b-', label='$h_\\times$')
  plt.legend()
  
  plt.figure('EOBNRv2 waveform in FD')
  plt.title('EOBNRv2 waveform in FD')
  plt.plot(fff, np.abs(hptf), 'r-', label='$h_+$')
  plt.plot(fff, np.abs(hctf), 'b-', label='$h_\\times$')
  plt.legend()
  
  plt.figure('TD waveform and noise')
  plt.title('TD waveform and noise')
  plt.plot(tt, hpt, 'r-', label='EOBNRv2 h_+')
  plt.plot(tt[:len(noise_TD)], noise_TD, 'b.', label='AdvLIGO noise')
  plt.plot(tt[:len(noise_TD)], hpt[:len(noise_TD)]+noise_TD, 'g.', label='Simulated signal')
  plt.legend()
  
  plt.figure('IMRPhenomB waveform in FD')
  plt.title('IMRPhenomB waveform in FD')
  plt.plot(ff, np.abs(hpf), 'r-', label='$h_+$')
  plt.plot(ff, np.abs(hcf), 'b-', label='$h_\\times$')
  plt.legend()
  
  #plt.figure('IMRPhenomB waveform in TD')
  #plt.title('IMRPhenomB waveform in TD')
  #plt.plot(np.fft.irfft(hpf), 'r-', label='$h_+$')
  #plt.plot(np.fft.irfft(hcf), 'b-', label='$h_\\times$')
  #plt.legend()

  plt.show()
