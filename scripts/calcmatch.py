""" 
Calculate the matches of NRAR waveforms with analytical templates 

P. Ajith, 2013-02-01 

$Id: calcmatch.py 165 2013-10-11 10:55:21Z ajith_p $
"""

import lal
import lalsimulation as lalsim
import lalinspiral as lalinsp 
from lal import MSUN_SI, MTSUN_SI, PC_SI, PI, PC_SI, C_SI 
import lalinspiral.sbank.waveforms as wf
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as P
import numpy as np 
import overlaps as ovp
import wavegen as wg 
import os, socket  
from datetime import datetime
import chirptimes as ct 
import imrtgrutils as tgr
import nr_fits as nr


##  input parameters 
resultsDirRoot = 'IMRTGR_WaveformSystStudy_2017-09-11'				# directory to which results should be saved 
sampl_rate = 4096. 
psd = 'aLIGOZeroDetHighPower'
run_type = 'inspiral' 
#run_type = 'post-inspiral' 
N_templ = 1000 

# parameters of the target waveform
m1_targ = 24. 
m2_targ = 25. 
spin1_targ  = [0., 0., 0.]
spin2_targ  = [0., 0., 0.]
domain_targ = 'time'
approx_targ = 'SEOBNRv2'
phaseO_targ = -1 
amplO_targ = -1 

q_targ = 1.15 
m_targ = 60. 
targ_wave_loc = '/Users/pajith/working/cbc/testGR_IR/data/IHES_EOB_mod/waveforms/popsynth_modGR'
targ_wave_file_name = '%s/injection_1013_q_1.15_a1_1_.dat' %targ_wave_loc
approx_targ = 'NRPNHybridtxt'
m1_targ = m_targ/(1.+ q_targ) 
m2_targ = q_targ*m_targ/(1.+q_targ) 
dt_M = 0.1750
sampl_rate = 1./(dt_M*m_targ*lal.MTSUN_SI)

Mf_targ, chif_targ = tgr.calc_final_mass_spin(m1_targ, m2_targ, 0., 0., 0, 0, 0., 0., 0., fit_formula='nonprecspin_Healy2014')
f_isco_Kerr = nr.calc_isco_freq(chif_targ)/(Mf_targ*lal.MTSUN_SI)

if run_type == 'inspiral': 
	f_min = 20. 															# low freq cutoff to be used in the match calculation 
	f_max = f_isco_Kerr															# high freq cutoft to be used in the match calculation 
elif run_type == 'post-inspiral': 
	f_min = f_isco_Kerr															# high freq cutoft to be used in the match calculation 
	f_max = sampl_rate/2. 

# parameters of the template waveform 
m1_templ = 24.
m2_templ = 25. 
spin1_templ = [0., 0., 0.]
spin2_templ = [0., 0., 0.]
domain_templ = 'freq'
approx_templ = 'SEOBNRv2_ROM_DoubleSpin'
#approx_templ = 'IMRPhenomPv2'
amplO_templ = -1 
phaseO_templ = -1 

# sky location and orientation of the binary in the detector frame. For optimally 
# located and oriented binary, set everything to zero 
sky_theta = 0
sky_phi = 0
polarization_angl = 0
incl_angle = 0
tag = '%s_vs_%s_%s' %(approx_targ, approx_templ, run_type) 

plot_fft = 0
print_waveforms = 0

# copy the scipt to the output directory 
resultsDirRoot = '%s/%s' %(resultsDirRoot, tag)
os.system('mkdir -p %s' %resultsDirRoot)
os.system('cp %s %s/' %(__file__, resultsDirRoot))

# create results directory
resultsDir = '%s/%s' %(resultsDirRoot, psd)
os.system('mkdir -p %s' %(resultsDir))

fp = open('%s/Match_%s_%s_m1_%.2f_m2_%.2f_spintempl_extended_eta_2.dat' %(resultsDir, tag, psd, m1_targ, m2_targ), 'a+')
fp.write('# Produced by: %s\n' %os.path.basename(__file__))
fp.write('# Date: %s\n' %datetime.now().isoformat()) 
fp.write('# Host: %s\n' %socket.gethostname())
fp.write('# m1  m2  s1x  s1y  s1z  s2x  s2y  s2z  f_min  f_max  match  \n')

# we distinguish between the lower and upper frequency cutoff used to generate the 
# waveform (f_min_wave and f_max_wave) from the lower and upper frequency cutoffs 
# (f_min and f_max) used in the integraation to compute the match function 
f_min_wave = min([f_min/2., 20.]) 
f_max_wave = sampl_rate/2. 

# crate a range of template masses 
m_targ = m1_targ+m2_targ 
q_targ = m1_targ/m2_targ

m_templ_vec = np.random.uniform(m_targ*0.9, m_targ*1.1, N_templ)
eta_templ_vec = np.random.uniform(0.20, 0.221, N_templ)
spin1z_templ_vec = np.random.uniform(-0.98, 0.98, N_templ)
spin2z_templ_vec = np.random.uniform(-0.98, 0.98, N_templ)

spin1x, spin1y, spin1z, spin2x, spin2y, spin2z = 0., 0., 0., 0., 0., 0.



for i_mass in range(len(m_templ_vec)):

	m_templ = m_templ_vec[i_mass]
	eta_templ = eta_templ_vec[i_mass] 
	spin1z = spin1z_templ = spin1z_templ_vec[i_mass]
	spin2z = spin2z_templ = spin2z_templ_vec[i_mass]

	m1_templ, m2_templ = ovp.individual_masses(m_templ, eta_templ)

	# compute the chirp time and the number of samples required 
	m1 = min([m1_targ, m1_templ])
	m2 = min([m2_targ, m2_templ])
	#tau = ct.calc_chirptime(m1, m2, f_min_wave)
	#N = int(2**np.ceil(np.log2(tau*sampl_rate)))
	N = 599269

	# calculate hte match  
	# match, snr_sq = ovp.calc_match_waveforms(m1_targ, m2_targ, spin1_targ, spin2_targ, m1_templ, m2_templ, spin1_templ, spin2_templ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min, f_max, f_min_wave, f_max_wave, sampl_rate, N, approx_targ, amplO_targ, phaseO_targ, approx_templ, amplO_templ, phaseO_templ, domain_targ, domain_templ, psd, plot_fft, print_waveforms, targ_wave_file_name=0)

	# calculate match - OLD CALCFF CODE 
	match = ovp.calc_match_waveforms(m1_targ, m2_targ, m1_templ, m2_templ, spin1x, spin1y, spin1z_templ, spin2x, spin2y, spin2z_templ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min, f_max, f_min_wave, f_max_wave, sampl_rate, N, approx_targ, amplO_targ, phaseO_targ, approx_templ, amplO_templ, phaseO_templ, domain_targ, domain_templ, psd, plot_fft, print_waveforms, targ_wave_file_name)

	# save the plot of the targ and templ waveform 
	if plot_fft: 
		os.system('mv fft_whitened_*.png %s' %resultsDir)

	# save the restuls 
	#spin1x, spin1y, spin1z = spin1_templ
	#spin2x, spin2y, spin2z = spin2_templ
	fp.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.6e\n' %(m1_templ, m2_templ, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_min, f_max, match))
	print '... m1 = %2.1f m2 = %2.1f s1 = [%3.2f, %3.2f, %3.2f] s2 = [%3.2f, %3.2f, %3.2f] f = [%3.2f, %3.2f] match = %5.4f' %(m1_templ, m2_templ, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_min, f_max, match)

fp.close() 

