"""
Code to generate GR waveform
"""

from numpy import sqrt, sin, cos, pi,exp
import matplotlib
matplotlib.use('Agg')
import sys, numpy as np
#import matplotlib.pyplot as plt
import phenomhh as phh
from lal import MSUN_SI, MTSUN_SI, PC_SI, PI, PC_SI, C_SI, GAMMA, MRSUN_SI
import pycbc.filter.matchedfilter as mfilter
import pycbc.psd
import pycbc.noise.gaussian
from pycbc  import  detector
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pycbc.types.frequencyseries
from pycbc.types import TimeSeries, FrequencySeries, zeros
sys.path.insert(0, '/home/abhirup/Documents/Work/testGR_IR/src')
import template_lalsiminspiral as phhsi_lal
import standard_gwtransf as gw

###### MAIN ##########
### parameters of the output data####
f_low=20.
dt = 0.0005

m1, m2 = 36., 29.
M = m1 + m2          ###### Total Mass in M_SUN
Mc, q = gw.mcq_from_comp(m1, m2)
r = 500.
iota=0.
Psi_ref=0.
t0=0.          ###### time of arrival

ra=0.          ##### Sky localisation
dec =0.
pol=0.

approx = 'IMRPhenomPv2'
################################################


t, hpt, hct = phhsi_lal.lalsiminspiraltd_waveform_SI(approx, Mc, q, r, iota, Psi_ref, f_low, dt)
df = 1./(len(hpt)*dt)

Fp,Fc = detector.overhead_antenna_pattern(ra, dec, pol)
psd = pycbc.psd.aLIGOZeroDetHighPower(int((len(hpt)+1)/2), df, f_low)

signal=Fp*hpt+Fc*hct
signal_time=pycbc.types.timeseries.TimeSeries(signal,delta_t=dt)
SNR=mfilter.sigma(signal_time,psd=psd,low_frequency_cutoff=f_low)

print 'The SNR is... %f'%SNR
print 'dL for the output SNR is.. %f Mpc'%r
print 'Chirp mass is.. %f solar mass'%Mc
print 'Lower cutoff.. %f Hz'%f_low

##### Output location and data file name #######

'''
Description of output data:

Output data contains [freq ,Real.part of signal ,Imag.part of signal ,PSD] in the same order
Data is non-zero starting from f_low. Last entry of data is (N-2)*df

'''
out_file = 'GR_M_%.2f_Mc_%.2f_q_%.2f_snr_%.2f_i_%.2f_dL_%.2f_phiref_%.2f_t0_%.2f_ra_%.2f_dec_%.2f_psi_%.2f_flow_20Hz_dt_%f_TD.dat'%(M,Mc,q,SNR,iota, r, Psi_ref, t0, ra, dec, pol, dt)

### Generating noise from psd ###
#noise=pycbc.noise.gaussian.frequency_noise_from_psd(psd, seed=None)
data=signal#+noise ## comment noise to generate noise free data
#################################

np.savetxt(out_file, np.c_[t,data, hpt, hct],header='t data hpt hct')
