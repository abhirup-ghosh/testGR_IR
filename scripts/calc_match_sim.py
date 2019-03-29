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
from matplotlib import pyplot as P
import numpy as np 
import overlaps as ovp
import wavegen as wg 
import os, socket  
from datetime import datetime
from optparse import OptionParser



##  input parameters 
f_min = 30. 								# low freq cutoff of the detector (Hz) 
resultsDirRoot = 'Test-2015-10-12'				# directory to which results should be saved 
tag = 'SEOBNRv2_ROM_DS_vs_SEOBNRv2_ROM_DS'
overlapType = 'max'	

os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
psd_path = os.path.join(home, 'src/lalsuite/lalsimulation/src')							# overlap type is max
psd = os.path.join(psd_path, 'LIGO-P1200087-v18-aLIGO_EARLY_HIGH.txt')					# detector noise 
m1 = 10											# list of total masses 
m2 = 40											# list of total masses 
chi = 0.
sampl_rate = 2048. 
N = sampl_rate*32. ###FIXME
spin1x = spin1y = spin2x = spin2y = 0 
det_name = 'aLIGO_EARLY_HIGH'


domain_targ = 'freq'
domain_templ = 'freq'
approx_targ = 'SEOBNRv2_ROM_DoubleSpinthreePointFivePN'
approx_templ = 'SEOBNRv2_ROM_DoubleSpinthreePointFivePN'
amplO_templ = amplO_targ = -1
phaseO_templ = phaseO_targ = 7 ###FIXME

# sky location and orientation of the binary in the detector frame. For optimally 
# located and oriented binary, set everything to zero 
sky_theta = 0
sky_phi = 0
polarization_angl = 0
incl_angle = 0

plot_fft = 0
print_waveforms = 0

# copy the scipt to the output directory 
#resultsDirRoot = '%s/%s' %(resultsDirRoot, tag)
#os.system('mkdir -p %s' %resultsDirRoot)
#os.system('cp %s %s/' %(__file__, resultsDirRoot))

# create results directory
#resultsDir = '%s/chi%2.1f/%s' %(resultsDirRoot, chi, psd)
#os.system('mkdir -p %s' %(resultsDir))
#outFileTag = '%s/%s' %(resultsDir, tag)

#fp = open('%s_FFPhysParams_%s.dat' %(outFileTag, psd), 'w+')
#fp.write('# Produced by: %s\n' %os.path.basename(__file__))
#fp.write('# Date: %s\n' %datetime.now().isoformat()) 
#fp.write('# Host: %s\n' %socket.gethostname())
#fp.write('# m1  m2  s1x  s1y  s1z  s2x  s2y  s2z  f_min  f_max  match  \n')

m = m1+m2 
eta = m1*m2/m**2.
spin1z = spin2z = chi/(1.-76.*eta/113.)
f_max = (1./6.)**(3./2)/(np.pi*m*MTSUN_SI)

# we distinguish between the lower and upper frequency cutoff used to generate the 
# waveform (f_min_wave and f_max_wave) from the lower and upper frequency cutoffs 
# (f_min and f_max) used in the integraation to compute the match function 
f_min_wave = f_min
f_max_wave = sampl_rate/2. 

m1_targ = m1
m2_targ = m2

N = 100
m1_templ_list = np.random.uniform(1., 100., N)
m2_templ_list = np.random.uniform(1., 100., N)

spin1_targ = [spin1x, spin1y, spin1z]
spin2_targ = [spin2x, spin2y, spin2z]
spin1_templ = [spin1x, spin1y, spin1z]
spin2_templ = [spin2x, spin2y, spin2z]

match = np.zeros_like(m1_templ_list)

for i in range(N):

	m1_templ = m1_templ_list[i]
	m2_templ = m1_templ_list[i]
	
	match[i] = ovp.calc_match_waveforms(m1_targ, m2_targ, spin1_targ, spin2_targ, m1_templ, m2_templ, spin1_templ, spin2_templ, incl_angle, sky_theta, sky_phi, polarization_angl, f_min, f_max, f_min_wave, f_max_wave, sampl_rate, N, approx_targ, amplO_targ, phaseO_targ, approx_templ, amplO_templ, phaseO_templ, domain_targ, domain_templ, det_name, plot_fft, print_waveforms, targ_wave_file_name=0)

P.figure()
P.scatter(m1_templ, m1_templ, s=20, c=match, mew=0)
P.colorbar()
P.xlabel('$m_1~[M_\odot]$')
P.xlabel('$m_2~[M_\odot]$')
P.plot(m1_targ, m2_targ, 'ko')
P.grid()
P.savefig('match.png')
P.show()

#if plot_fft: 
#	os.system('mv fft_whitened_targ_templ_%s.png %s/fft_targ_templ_%s_m%2.1f.png' %(tag, resultsDir, tag, m))

#fp.write('%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%5.4f\n' %(m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_min, f_max, match))
#print '... m1 = %2.1f m2 = %2.1f s1 = [%3.2f, %3.2f, %3.2f] s2 = [%3.2f, %3.2f, %3.2f] f = [%3.2f, %3.2f] match = %5.4f' %(m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_min, f_max, match)
print match
#fp.close() 
