""" 
Calculate the matches of NRAR waveforms with analytical templates 

P. Ajith, Abhirup Ghosh 2017-09-03 

$Id: calcmatch.py 165 2013-10-11 10:55:21Z ajith_p $
"""

import lal
import lalsimulation as lalsim
import lalinspiral as lalinsp 
from lal import MSUN_SI, MTSUN_SI, PC_SI, PI, PC_SI, C_SI 
import lalinspiral.sbank.waveforms as wf
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 
import overlaps as ovp
import wavegen as wg 
import os, socket  
from datetime import datetime
import chirptimes as ct 
import scipy
from scipy import interpolate
from optparse import OptionParser
import imrtestgr as tgr
import nr_fits as nr

parser = OptionParser()
parser.add_option("--m1-targ", dest="m1_targ", type='float')
parser.add_option("--m2-targ", dest="m2_targ", type='float')
parser.add_option("--portion", dest="portion", type='string')
parser.add_option("--approx-targ", dest="approx_targ", type='string')
parser.add_option("--domain-targ", dest="domain_targ", type='string')
parser.add_option("--approx-templ", dest="approx_templ", type='string')
parser.add_option("--domain-templ", dest="domain_templ", type='string')
(options, args) = parser.parse_args()

m1_targ = options.m1_targ
m2_targ = options.m2_targ
portion = options.portion
approx_targ = options.approx_targ
domain_targ = options.domain_targ
approx_templ = options.approx_templ
domain_templ = options.domain_templ

Mf, af = tgr.calc_final_mass_spin(m1_targ, m2_targ, 0., 0., 'nonprecspin_Healy2014')
f_isco = nr.calc_isco_freq(af)/(Mf*MTSUN_SI)

if portion == 'IMR':
	  f_min, f_max = 20., 2048.
if portion == 'inspiral':
          f_min, f_max = 20., f_isco
if portion == 'post-inspiral':
          f_min, f_max = f_isco, 2048.

##  input parameters 
resultsDirRoot = '../runs/waveform_systematics/IMRTGR_WaveformSystStudy_2017-08-31_m1m2_0.99-1.01_a12z_-1.0-1.0' # directory to which results should be saved 
sampl_rate = 4096. 
psd = 'aLIGOZeroDetHighPower'

# parameters of the target waveform
spin1_targ  = [0., 0., 0.]
spin2_targ  = [0., 0., 0.]
phaseO_targ = -1 
amplO_targ = -1 

# parameters of the template waveform 
amplO_templ = -1 
phaseO_templ = -1 

# sky location and orientation of the binary in the detector frame. For optimally 
# located and oriented binary, set everything to zero 
sky_theta = 0
sky_phi = 0
polarization_angl = 0
incl_angle = 0
tag = '%s_vs_%s' %(approx_targ, approx_templ) 

plot_fft = 0
print_waveforms = 0

# copy the scipt to the output directory 
resultsDirRoot = '%s/%s' %(resultsDirRoot, tag)
os.system('mkdir -p %s' %resultsDirRoot)
os.system('cp %s %s/' %(__file__, resultsDirRoot))

# create results directory
resultsDir = '%s/%s/%s' %(resultsDirRoot, psd, portion)
os.system('mkdir -p %s' %(resultsDir))

fp = open('%s/Match_%s_%s_m1_%.2f_m2_%.2f_a1z_%.2f_a2z_%.2f.dat' %(resultsDir, tag, psd, m1_targ, m2_targ, spin1_targ[2], spin2_targ[2]), 'a+')
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
m1_templ_vec = np.random.uniform(0.99*m1_targ, 1.01*m1_targ, 1000)
m2_templ_vec = np.random.uniform(0.99*m2_targ, 1.01*m2_targ, 1000)
spin1z_templ_vec = np.random.uniform(-1., 0.99, 1000)
spin2z_templ_vec = np.random.uniform(-1., 0.99, 1000)

for i_mass in range(len(m1_templ_vec)):

	m1_templ = m1_templ_vec[i_mass]
	m2_templ = m2_templ_vec[i_mass]
	spin1_templ = [0., 0., spin1z_templ_vec[i_mass]]
	spin2_templ = [0., 0., spin2z_templ_vec[i_mass]]

	# compute the chirp time and the number of samples required 
	m1 = min([m1_targ, m1_templ])
	m2 = min([m2_targ, m2_templ])
	tau = ct.calc_chirptime(m1, m2, f_min_wave)
	N = int(2**np.ceil(np.log2(tau*sampl_rate)))

	# calculate hte match 
	match, snr_sq = ovp.calc_match_waveforms(m1_targ, m2_targ, spin1_targ, spin2_targ, m1_templ, m2_templ, spin1_templ, spin2_templ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min, f_max, f_min_wave, f_max_wave, sampl_rate, N, approx_targ, amplO_targ, phaseO_targ, approx_templ, amplO_templ, phaseO_templ, domain_targ, domain_templ, psd, plot_fft, print_waveforms, targ_wave_file_name=0)

	# save the plot of the targ and templ waveform 
	if plot_fft: 
		os.system('mv fft_whitened_*.png %s' %resultsDir)

	# save the restuls 
	spin1x, spin1y, spin1z = spin1_templ
	spin2x, spin2y, spin2z = spin2_templ
	fp.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.6e\n' %(m1_templ, m2_templ, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_min, f_max, match))

fp.close() 

"""
data = np.genfromtxt('%s/Match_%s_%s_m1_%.2f_m2_%.2f.dat' %(resultsDir, tag, psd, m1_targ, m2_targ), names=True, dtype=None, skip_header=3)
m1 = data['m1']
m2 = data['m2']
match = data['match']
m1i = np.linspace(m1.min(), m1.max(), len(m1))
m2i = np.linspace(m2.min(), m2.max(), len(m2))
matchi = scipy.interpolate.griddata((m1, m2), match, (m1i[None,:], m2i[:,None]), method='cubic')

plt.figure()
plt.scatter(m1, m2, c=np.log10(1.-match), lw=0, alpha=1)
plt.colorbar()
CS = plt.contour(m1i, m2i, matchi, levels=(0.68, 0.95))
plt.clabel(CS, inline=1, fontsize=10)
plt.axvline(x=m1_targ)
plt.axhline(y=m2_targ)
plt.xlabel('$m_1$')
plt.ylabel('$m_2$')
plt.title('$m_1$=%.2f; $m_2$=%.2f; colorbar: log$_{10}$(match)'%(m1_targ, m2_targ))
plt.grid()
plt.xlim([0.99*m1_targ, 1.01*m1_targ])
plt.ylim([0.99*m2_targ, 1.01*m2_targ])
plt.savefig('%s/Match_%s_%s_m1_%.2f_m2_%.2f.png' %(resultsDir, tag, psd, m1_targ, m2_targ))
"""
