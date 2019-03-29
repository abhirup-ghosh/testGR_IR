#!/usr/bin/env python

"""
This script creates an injection and a .ini file to run LALInference

(C) Archisman Ghosh, Abhirup Ghosh, 2015-09-20
"""

import sys, os, commands, numpy as np
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
import imrtestgr as tgr
import nr_fits as nr
import matplotlib.pyplot as plt

os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
lal_prefix = os.getenv('LAL_PREFIX')
psd_path = '%s/share/lalsimulation'%lal_prefix
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')

# ### INJECTION ###
no_inj = 10000

# Waveform info
approximant = 'SEOBNRv4pseudoFourPN' # injection waveform 
amp_order = 0 # amplitude PN order of the injection waveform
f_low = 8. # low-frequency cutoff (Hz)
waveform_info = '--waveform %s --amp-order %d --f-lower %f --taper-injection startend'%(approximant, amp_order, f_low)

# PSD info
f_start = 10.
ligo_psd = 'LALAdLIGO'
virgo_psd = 'LALAdVirgo'
psd_info = ' '#'--ifos H1,L1,V1 --ligo-fake-psd %s --ligo-start-freq %f --virgo-fake-psd %s --virgo-start-freq %f'%(ligo_psd, f_start, virgo_psd, f_start)

# Time info
gps_start_time = 1126285216 # O1 start time 
time_step = 2630./np.pi # time interval between nearby injections 
gps_end_time = gps_start_time + time_step*no_inj
time_info = '--t-distr fixed --gps-start-time %f --gps-end-time %f --time-step %f'%(gps_start_time, gps_end_time, time_step)

# Parameters
m_min = 10.
m_max = 80.
a_min = 0.
a_max = 1.
phi_c = 0.
psi = 0.
z_min = 0.
z_max = 0.5

param_info = '--m-distr componentMass --min-mass1 %f --max-mass1 %f --min-mass2 %f --max-mass2 %f --enable-spin --aligned --min-spin1 %f --max-spin1 %f --min-spin2 %f --max-spin2 %f --i-distr uniform --coa-phase-distr fixed --fixed-coa-phase %f --polarization %f --l-distr random --d-distr volume --min-z %f --max-z %f'%(m_min, m_max, m_min, m_max, a_min, a_max, a_min, a_max, phi_c, psi, z_min, z_max)

# Output info
inj_file = '../injections/waveform_systematics/waveform_systematics_investigations_%s_m_10_80_z_0_0p5'%(approximant) # output file name for xml and txt file
output_info = '--output %s.xml'%(inj_file)

# Generate the injection xml file by calling lalapps_inspinj 
run_command = 'lalapps_inspinj %s %s %s %s %s'%(waveform_info, psd_info, time_info, param_info, output_info)
print run_command
os.system('echo %s > %s.out'%(run_command, inj_file))
os.system(run_command)

os.system('ligolw_print -d " " -c mass1 -c mass2 -c spin1x -c spin2x -c spin1y -c spin2y -c spin1z -c spin2z -c longitude -c latitude -c inclination -c polarization -c distance -c geocent_end_time -c geocent_end_time_ns -c coa_phase %s.xml > %s_all_param.tmp' %(inj_file, inj_file))
os.system('echo \# mass1 mass2 spin1x spin2x spin1y spin2y spin1z spin2z longitude latitude inclination polarization distance geocent_end_time geocent_end_time_ns coa_phase > header.txt')
os.system('cat header.txt %s_all_param.tmp > %s.txt' %(inj_file, inj_file))
os.system('rm %s_all_param.tmp header.txt' %(inj_file))
