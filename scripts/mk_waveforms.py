""" 
Generate waveforms with transparency weighted by the logL values from posterior samples 

P. Ajith, 2015-11-11 

$Id:$
"""

import lalsimulation as lalsim
from lal import MSUN_SI, MTSUN_SI, PC_SI, PI, PC_SI, C_SI 
from matplotlib import pyplot as P
import numpy as np, os
from datetime import datetime
import time  
import plotsettings 
import imrtestgr as tgr
import nr_fits as nr



##  input parameters 
fs = 4096. 								# sampling rate of the data 
spin1x = 0.								# spin components 
spin1y = 0.
spin1z = 0.
spin2x = 0.
spin2y = 0 
spin2z = 0.
approx = 'SEOBNRv2'					# name of the PN approximant 
amplO = 0							# PN order of the amplitude (7 = 3.5PN) 
phaseO = 8 							# PN order of the phase 
f_low = 30.
incl_angle = 0.						# inclination angle 
f_final = 0.
t_max = 0.22	
spin_fit_formula = 'nonprecspin_Healy2014'

t_ref = 0.
m1_vec = [25., 25., 25., 25., 25]
m2_vec = [25., 5., 25., 25., 25.]
spin1z_vec = [0., 0., -0.98, 0, 0.98]
spin2z_vec = [0., 0., -0.98, 0, 0.98]
	
# get the approximant enum for LALSimulation 
approx_enum = lalsim.GetApproximantFromString(approx)

#for i in range(len(m1_vec)):
for i in range(len(m1_vec)):

	# compute component masses, symm. mass ratio, f-lower etc 
	m1 = m1_vec[i]
	m2 = m2_vec[i]
	spin1z = spin1z_vec[i]
	spin2z = spin2z_vec[i]
	logl = 1.
	alph = 1.
	t0 = 0.
	dist = 1.
	m = m1+m2 
	q = m1/m2

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

	t = np.linspace(0, hp.data.length/fs, hp.data.length)+t0-t_ref

	end_time=time.time()


	# calculate the mass and spin of the final BH 
	Mf, af = tgr.calc_final_mass_spin(m1, m2, spin1z, spin2z, spin_fit_formula)

	# calculate the Schwarzschild ISCO 
	f_isco_Schw = 6**-1.5/(np.pi*(m1+m2)*MTSUN_SI) 

	# calculate the Kerr ISCO freq 
	f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*MTSUN_SI)

	# calculate the dominant QNM freq 
	f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*MTSUN_SI)
	
	# compute F(t) and find the index corresponding to F(f) = 134 
	hh = hp.data.data+1j*hc.data.data 
	Foft = np.gradient(np.unwrap(np.angle(hh)))/np.gradient(t)/(2*np.pi)
	f_isco_kerr_idx = np.where(Foft >= f_isco_Kerr)[0][0]
	f_isco_schw_idx = np.where(Foft >= f_isco_Schw)[0][0]
	#f_qnm_idx = np.where(Foft >= f_qnm)[0][0]
		

	P.figure(figsize=(8,4))
	P.subplot2grid((2,3), (0,0), colspan=2)
	P.plot(t, hp.data.data, color='k', alpha=alph, lw=1.5)
	P.plot(t, np.sqrt(hp.data.data**2+hc.data.data**2), color='orange', alpha=alph, lw=0.5)
	P.axvline(x=t[f_isco_kerr_idx], color='r', lw=0.5, ls='-')
	P.axvline(x=t[f_isco_schw_idx], color='m', lw=0.5, ls='-')
	P.ylabel('$h_+(t)$')
	P.title('$M = %2.1f M_\odot ~ q = %2.1f ~ \chi_1 = %3.2f ~ \chi_2 = %3.2f$' %(m, q, spin1z, spin2z), fontsize=10)

	P.subplot2grid((2,3), (0,2), colspan=1)
	P.plot(t, hp.data.data, color='k', alpha=alph, lw=1.5)
	P.plot(t, np.sqrt(hp.data.data**2+hc.data.data**2), color='orange', alpha=alph, lw=0.5)
	P.axvline(x=t[f_isco_kerr_idx], color='r', lw=0.5, ls='-')
	P.axvline(x=t[f_isco_schw_idx], color='m', lw=0.5, ls='-')
	P.xlim(t[f_isco_schw_idx]-0.01, max(t))

	P.subplot2grid((2,3), (1,0), colspan=2)
	P.plot(t, Foft, color='k', alpha=alph, lw=1.5)
	P.axvline(x=t[f_isco_kerr_idx], color='r', lw=0.5, ls='-')
	P.axvline(x=t[f_isco_schw_idx], color='m', lw=0.5, ls='-')
	P.ylabel('$F(t)$ [Hz]')
	P.xlabel('$t$ [sec]')
	P.ylim(20, 1.1*f_qnm)
	P.axhline(y=f_isco_Schw, label='Schw ISCO', color='m', lw=0.5, ls='-')
	P.axhline(y=f_isco_Kerr, label='Kerr ISCO', color='r', lw=0.5, ls='-')
	P.axhline(y=f_qnm, label='QNM', color='c', lw=0.5, ls='-')
	P.text(0.01, f_isco_Kerr+10, 'Kerr ISCO', fontsize=8, color='r')
	P.text(0.01, f_isco_Schw+10, 'Schw ISCO', fontsize=8, color='m')
	P.text(0.01, f_qnm+10, 'QNM', fontsize=8, color='c')

	P.subplot2grid((2,3), (1,2), colspan=2)
	P.plot(t, Foft, color='k', alpha=alph, lw=1.5)
	P.axvline(x=t[f_isco_kerr_idx], color='r', lw=0.5, ls='-')
	P.axvline(x=t[f_isco_schw_idx], color='m', lw=0.5, ls='-')
	P.xlabel('$t$ [sec]')
	P.ylim(20, 1.1*f_qnm)
	P.axhline(y=f_isco_Schw, label='Schw ISCO', color='m', lw=0.5, ls='-')
	P.axhline(y=f_isco_Kerr, label='Kerr ISCO', color='r', lw=0.5, ls='-')
	P.axhline(y=f_qnm, label='QNM', color='c', lw=0.5, ls='-')
	P.xlim(t[f_isco_schw_idx]-0.01, max(t))


	P.savefig('waveform_m%2.1f_q%2.1f_chi1%2.1f_chi2%2.1f.png' %(m, q, spin1z, spin2z), dpi=200)




