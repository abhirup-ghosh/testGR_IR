import matplotlib as mpl
mpl.use('Agg')
import os, lal, subprocess, scipy, numpy as np
from scipy.interpolate import interp1d
import standard_gwtransf as gw
import pycbc
from pycbc import waveform
import os
import matplotlib.pyplot as plt

def lalsiminspiral_waveform_SI(approximant,mc,q,dL,iota,t0,phi_c,flow,df,f_high):

  try:
	m1, m2 = gw.comp_from_mcq(mc, q)
	eta = gw.eta_from_q(q)
	mtot = m1 + m2
	hpf, hcf = waveform.get_fd_waveform(approximant=approximant,mass1=m1, mass2=m2, distance=dL, inclination=iota, coa_phase=phi_c, delta_f=df, f_lower=flow, f_final=f_high)

	print '... generated waveform for m1=%.2f, m2=%.2f, Mc=%.2f, q=%.2f, mtot=%.2f, eta=%.4f'%(m1,m2,mc,q, mtot, eta)

	f = np.array(hpf.sample_frequencies)
	hpf = (np.array(hpf.data))*np.exp(-2*np.pi*1j*f*t0)
	hcf = (np.array(hcf.data))*np.exp(-2*np.pi*1j*f*t0)

	return f, hpf, hcf
  except RuntimeError:
	print '... not generated waveform for m1=%.2f and m2=%.2f'%(m1,m2)
  except TypeError:
	print '... not generated waveform for m1=%.2f and m2=%.2f'%(m1,m2)



def lalsiminspiraltd_waveform_SI(approximant,mc,q,dL,iota,phi_c,flow,dt):

  try:
        m1, m2 = gw.comp_from_mcq(mc, q)
        eta = gw.eta_from_q(q)
        mtot = m1 + m2
        hpt, hct = waveform.get_td_waveform(approximant=approximant,mass1=m1, mass2=m2, distance=dL, inclination=iota, coa_phase=phi_c, delta_t=dt, f_lower=flow)

        print '... generated waveform for m1=%.2f, m2=%.2f, Mc=%.2f, q=%.2f, mtot=%.2f, eta=%.4f'%(m1,m2,mc,q, mtot, eta)

        t = hpt.sample_times
	t = np.array(t)
	hpt = np.array(hpt.data)
	hct = np.array(hct.data)

        return t, hpt, hct
  except RuntimeError:
        print '... not generated waveform for m1=%.2f and m2=%.2f'%(m1,m2)
  except TypeError:
        print '... not generated waveform for m1=%.2f and m2=%.2f'%(m1,m2)
