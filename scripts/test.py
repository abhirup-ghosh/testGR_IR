import numpy as np
import matplotlib.pyplot as plt
import pycbc.waveform, pycbc.noise, pycbc.psd, pycbc.filter


hp, hc = pycbc.waveform.get_td_waveform(approximant='SEOBNRv4', mass1=36, mass2=29, distance=500., delta_t=1.0/4096,f_lower=20)
delta_t = hp.delta_t

delta_f = 1. / len(hp) / delta_t

psd = pycbc.psd.analytical.aLIGOQuantumZeroDetHighPower(len(hp), delta_f, low_freq_cutoff=20.)

noise_td = pycbc.noise.gaussian.noise_from_psd(len(hp), delta_t, psd, seed=0)

snr = pycbc.filter.matchedfilter.sigma(hp, psd=psd, low_frequency_cutoff=20., high_frequency_cutoff=1024.)
print snr

