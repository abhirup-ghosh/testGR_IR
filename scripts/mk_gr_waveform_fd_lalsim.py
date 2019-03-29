"""
Code to generate GR waveform
"""

from numpy import sqrt, sin, cos, pi,exp
import matplotlib
matplotlib.use('Agg')
import numpy as np
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
import template_lalsiminspiral as phhsi_lal
import standard_gwtransf as gw

###### MAIN ##########
### parameters of the output data####
f_low=20.
df=0.1
N=10000

m1, m2 = 36, 29
M = m1 + m2          ###### Total Mass in M_SUN
Mc, q = gw.mcq_from_comp(m1, m2)
SNR_req=25.    ###### Required SNR
iota=pi/3.
Psi_ref=pi
t0=6.          ###### time of arrival

ra=1.          ##### Sky localisation
dec =1.
pol=0.

################################################


f, hpf, hcf = phhsi_lal.lalsiminspiral_waveform_SI(Mc, q, 1, iota, t0, Psi_ref, f_low, df, N)

Fp,Fc = detector.overhead_antenna_pattern(ra, dec, pol)
psd = pycbc.psd.aLIGOZeroDetHighPower(len(hpf), df, f_low)

signal=Fp*hpf+Fc*hcf
signal_freq=pycbc.types.frequencyseries.FrequencySeries(signal,delta_f=df,dtype=complex)
SNR=mfilter.sigma(signal_freq,psd=psd,low_frequency_cutoff=f_low)

print 'SNR at 1Mpc is... %f'%SNR 

r=SNR/SNR_req
f, hpf, hcf = phhsi_lal.lalsiminspiral_waveform_SI(Mc, q, r, iota, t0, Psi_ref, f_low, df, N)

signal=Fp*hpf+Fc*hcf
signal_freq=pycbc.types.frequencyseries.FrequencySeries(signal,delta_f=df,dtype=complex)
SNR=mfilter.sigma(signal_freq,psd=psd,low_frequency_cutoff=f_low)

Mc=(M*q**0.6)/((1.+q)**1.2)

print 'The SNR is... %f'%SNR
print 'dL for the output SNR is.. %f Mpc'%r
print 'Chirp mass is.. %f solar mass'%Mc

##### Output location and data file name #######

'''
Description of output data:

Output data contains [freq ,Real.part of signal ,Imag.part of signal ,PSD] in the same order
Data is non-zero starting from f_low. Last entry of data is (N-2)*df

'''
out_file = 'GR_M_%.2f_Mc_%.2f_q_%.2f_snr_%d_i_%.2f_dL_%.2f_phiref_%.2f_t0_%.2f_ra_%.2f_dec_%.2f_psi_%.2f_flow_20Hz.dat'%(M,Mc,q,iota,SNR, r, Psi_ref, t0, ra, dec, pol)

### Generating noise from psd ###
#noise=pycbc.noise.gaussian.frequency_noise_from_psd(psd, seed=None)
data=signal#+noise ## comment noise to generate noise free data
#################################

datar=np.real(data)
datai=np.imag(data)

np.savetxt(out_file, np.c_[f[:-1],datar[:-1],datai[:-1],psd[:-1]], header='f real_data imag_data psd')

