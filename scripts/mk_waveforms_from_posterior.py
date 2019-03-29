""" 
Generate waveforms with transparency weighted by the logL values from posterior samples 

P. Ajith, 2015-11-11 

$Id:$
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P
import lalsimulation as lalsim
from lal import MSUN_SI, MTSUN_SI, PC_SI, PI, PC_SI, C_SI 
import numpy as np, os
from datetime import datetime
import time  
import plotsettings 
from pylal import antenna


##  input parameters 
fs = 2048. 								# sampling rate of the data 
approx = 'SEOBNRv2'					# name of the PN approximant 
amplO = 0							# PN order of the amplitude (7 = 3.5PN) 
phaseO = 8 							# PN order of the phase 
f_low = 30.
f_final = 0.
t_max = 0.21	
spin1x = 0.
spin1y = 0.
spin2x = 0.
spin2y = 0.


post_file = '../data/GW150914/posterior_samples_imr.dat.gz'

# read posterior data 
data = np.genfromtxt(post_file, dtype=None, names=True)
m1_vec, m2_vec, logl_vec, tC_L_vec, tC_H_vec, ra_vec, dec_vec, dist_vec, incl_vec, psi_vec = data['m1'], data['m2'], data['logl'], data['l1_end_time'], data['h1_end_time'], data['ra'], data['dec'], data['dist'], data['theta_jn'], data['psi']
if ('a1' in data.dtype.names) and ('a2' in data.dtype.names):
  spin1z_vec, spin2z_vec = data['a1'], data['a2']
else:
  spin1z_vec, spin2z_vec = np.zeros(len(m1)), np.zeros(len(m2))

t_ref = np.min(tC_L_vec)

# sort all parameters according to the logL
idx = np.flipud(np.argsort(logl_vec))
logl_vec = logl_vec[idx]
m1_vec = m1_vec[idx]
m2_vec = m2_vec[idx]
spin1z_vec = spin1z_vec[idx]
spin2z_vec = spin2z_vec[idx]
tC_L_vec = tC_L_vec[idx]
tC_H_vec = tC_H_vec[idx]
dist_vec = dist_vec[idx]
ra_vec = ra_vec[idx]
dec_vec = dec_vec[idx]
incl_vec = incl_vec[idx]
psi_vec = psi_vec[idx]

# P.figure(figsize=(6,6))
# P.scatter(m1_vec, m2_vec, s=0.5, c=logl_vec, lw=0)
# P.colorbar()
# P.grid()
# P.savefig('posterior_samples_sorted.png', dpi=300)

P.figure(figsize=(12,8))

logl_max = np.max(logl_vec)
logl_min = np.min(logl_vec) 
	
alpha_vec = (logl_vec-logl_min)/(logl_max-logl_min)
alpha_vec = (np.exp(logl_vec-logl_max))**2
	
# get the approximant enum for LALSimulation 
approx_enum = lalsim.GetApproximantFromString(approx)

#for i in range(len(m1_vec)):
for i in range(10):

	# compute component masses, symm. mass ratio, f-lower etc 
	m1 = m1_vec[i]
	m2 = m2_vec[i]
	spin1z = spin1z_vec[i]
	spin2z = spin2z_vec[i]
	logl = logl_vec[i]
	alph = alpha_vec[i]
	tC_L = tC_L_vec[i]
	tC_H = tC_H_vec[i]
	dist = dist_vec[i]
	ra = ra_vec[i]
	dec = dec_vec[i]
	incl_angle = incl_vec[i]
	psi = psi_vec[i]

	start_time=time.time()

	# Generate template waveforms in the time domain. 
	hp, hc = lalsim.SimInspiralChooseTDWaveform(
		0.,						# phi_ref (rad) 
		1./fs,					# delta_t (sec) 
		m1 * MSUN_SI, 	# mass 1 (kg) 
		m2 * MSUN_SI,		# mass 2 (kg) 
		spin1x,					# spin1x (S1x/m1^2)  
		spin1y, 			# spin1y (S1y/m1^2)
		spin1z, 			# spin1z (S1z/m1^2)
		spin2x, 			# spin2x (S2x/m2^2)
		spin2y, 			# spin2y (S2y/m2^2)
		spin2z,					# spin2z (S2z/m2^2)
		f_low, 				# f_low (Hz) 
		f_low, 				# f_ref (Hz) 
		dist * 1.0e6 * PC_SI, 	# distance (m)
		incl_angle, 		# incl_angle angle (rad) 
		0.,						# lambda1 (Love number) 
		0.,						# lambda2 (Love number) 
		None, 				# waveform flags
		None, 				# non-GR parameters
		amplO, 				# amplitude PN order 
		phaseO, 			# phase PN order 
		approx_enum 		# approximant 
		)

	# unpack the data 
	N = hp.data.length
	hp = hp.data.data
	hc = hc.data.data

	# compute the antenna pattern functions
	Fp_H, Fc_H, Fav_H, Q_H = antenna.response(tC_H, ra, dec, incl_angle, psi, 'radians', 'H1' )
	Fp_L, Fc_L, Fav_L, Q_L = antenna.response(tC_L, ra, dec, incl_angle, psi, 'radians', 'H1' )

	# compute hoft at H1 and L1 
	h_H = Fp_H*hp + Fc_H*hc 
	h_L = Fp_L*hp + Fc_L*hc 
	t_L = np.linspace(0, N/fs, N)+tC_L-t_ref
	t_H = np.linspace(0, N/fs, N)+tC_H-t_ref

	# compute F(t) and find the index corresponding to F(f) = 134 
	hh = hp+1j*hc 
	#Foft = np.gradient(np.unwrap(np.angle(hh)))/np.gradient(t)/(2*np.pi)
	#f0idx = np.where(Foft >= 134.)[0][0]
	
	P.subplot(321)
	P.plot(t_H, h_H, color='grey', alpha=alph, lw=1)
	P.ylabel('$h(t)$')
	P.title('H1')
	P.subplot(322)
	P.plot(t_L, h_L, color='grey', alpha=alph, lw=1)
	P.ylabel('$h(t)$')
	P.title('L1')

# P.subplot(211)
# P.ylabel('$h(t)$')
# P.xlim(0, t_max)
# P.subplot(322)
# P.ylabel('$F(t)$ [Hz]')
# P.xlabel('$t$ [sec]')
# P.xlim(0, t_max)
# P.ylim(0, 300)
# P.axhline(y=134, color='r', lw=0.5, ls='-')
P.savefig('waveform.png', dpi=300)




