""" 
Compute the radiated energy as a function of the Fourier frequency (or time) from the waveforms 

P. Ajith, 2015-01-27
"""

import wavegen as wf 
import numpy as np 
import matplotlib.pyplot as plt 
import plotsettings 

def edf(m1, m2, f_low, fs, f_final):
	spin1x, spin1y, spin1z = 0., 0., 0.
	spin2x, spin2y, spin2z = 0., 0., 0.
	incl_angle, sky_theta, sky_phi, polarization_angl = 0., 0., 0., 0., 

#	f_low = 10. 
#	fs = 4096. 
#	f_final = fs/2.
	N = 16384*8
	approx = 'IMRPhenomB' #'EOBNRv2'
	dom='freq' # 'time'
	amplO, phaseO = 8, 8


	# generate the Fourier transform of h+(f) 
	hf = wf.generate_waveforms(m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, incl_angle, sky_theta, sky_phi, polarization_angl, f_low, f_final, fs, N, approx, amplO, phaseO, phi0=0., dataFileName=0, waveform_type=0, domain=dom, use_only_22=0, print_waveforms=0)

	# unpack the hf 
	hf = hf.data.data

	# construct a freq vector 
	f = np.linspace(0., f_final, len(hf)) 

	# compute the energy density and the cumulative energy loss 
	ef = (f*abs(hf))**2
	Ef = np.cumsum(ef)

	return f, hf, ef, Ef


f, hf, ef, Ef = edf(40., 10., 10, 2048., 1024.)

plt.figure()
plt.loglog(f, Ef/max(Ef))
plt.xlabel('$f$ [Hz]')
plt.grid()
plt.show()
