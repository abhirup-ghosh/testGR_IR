#!/usr/bin/python
""" 
Create the xml file of injections by running lalaapps_inspinj 

(C) Abhirup Ghosh, Archisman Ghosh & P. Ajith; Modified: 2014-09-27
$Id:$
"""
import os, numpy as np

# Executable
exec_script = 'lalapps_inspinj' # the lalapps_inspinj executable 

# Waveform info
#approximant = 'IMRPhenomBpseudoFourPN' #'EOBNRv2threePointFivePN' # injection waveform 
approximant = 'SEOBNRv2_ROM_DoubleSpinthreePointFivePN' # injection waveform 
amp_order = -1 # amplitude PN order of the injection waveform
f_low = 8 # low-frequency cutoff (Hz)
waveform_info = '--waveform %s --amp-order %d --f-lower %f --taper-injection startend'%(approximant, amp_order, f_low)

# PSD info
ligo_psd = 'LALAdLIGO' # psd to be assumed for ligo 
virgo_psd = 'LALAdVirgo' # psd to be assumed for virgo
psd_info = '' #'--ligo-fake-psd %s --virgo-fake-psd %s'%(ligo_psd, virgo_psd)

# Time info
gps_time = 1041033614 # gps time 
time_step = 2630./np.pi # time interval between nearby injections 
time_info = '--t-distr fixed --gps-start-time %f --gps-end-time %f --time-step %f'%(gps_time, gps_time, time_step)

# Parameters
M_list = [50., 75., 100., 150., 200.] # total mass (M_sun)
q_list = [1., 2., 4.]
#m1 = 15. # component mass 1 (M_sun)
#m2 = 15. # component mass 2 (M_sun)
#dL = 100*1000. # distance (kpc)
snr_min = 50
snr_max = snr_min
f_start = 30.
iota = 0. # inclination (degrees)
phi_c = 0. # coalescent phase (degrees)
psi = 0. # polarization (degrees)
alpha = 0. # right ascension (degrees)
delta = 0. # declination (degrees)
#param_info = '--m-distr fixMasses --fixed-mass1 %f --fixed-mass2 %f --d-distr uniform --min-distance %f --max-distance %f --i-distr fixed --fixed-inc %f --coa-phase-distr fixed --fixed-coa-phase %f --polarization %f --l-distr fixed --latitude %f --longitude %f --disable-spin'%(m1, m2, dL, dL, iota, phi_c, psi, alpha, delta)

for M in M_list:
  for q in q_list:
	m1 = M/(1.+q)
	m2 = M*q/(1.+q)
	print m1, m2
	#param_info = '--m-distr fixMasses --fixed-mass1 %f --fixed-mass2 %f --d-distr uniform --min-distance %f --max-distance %f --i-distr fixed --fixed-inc %f --coa-phase-distr fixed --fixed-coa-phase %f --polarization %f --l-distr fixed --latitude %f --longitude %f --disable-spin'%(m1, m2, dL, dL, iota, phi_c, psi, alpha, delta)
	param_info = '--m-distr fixMasses --fixed-mass1 %f --fixed-mass2 %f --snr-distr uniform --min-snr %f --max-snr %f --ifos H1 --ligo-start-freq %f --ligo-fake-psd %s --i-distr fixed --fixed-inc %f --coa-phase-distr fixed --fixed-coa-phase %f --polarization %f --l-distr fixed --latitude %f --longitude %f --disable-spin'%(m1, m2, snr_min, snr_max, f_start, ligo_psd, iota, phi_c, psi, alpha, delta)

# Output info
	#out_file = '../injections/example_inj_IMRPhenomB_snr50_%s_%s'%(M,q) # output file name for xml and txt file
	out_file = '../injections/example_inj_SEOBNRv2_ROM_DS_snr50_%s_%s'%(M,q) # output file name for xml and txt file
	output_info = '--output %s.xml'%(out_file)

# Generate the injection xml file by calling lalapps_inspinj 
	run_command = '%s %s %s %s %s %s'%(exec_script, waveform_info, psd_info, time_info, param_info, output_info)
	print(run_command)
	os.system(run_command)

# Print some relevant columns of the xml file to an ASCII table with header 
	os.system('ligolw_print -d " " -c mass1 -c mass2 -c longitude -c latitude -c inclination -c polarization -c distance -c geocent_end_time -c geocent_end_time_ns -c coa_phase %s.xml > %s.tmp' %(out_file, out_file))
	os.system('echo \# mass1 mass2 longitude latitude inclination polarization distance geocent_end_time geocent_end_time_ns coa_phase > header.txt')
	os.system('cat header.txt %s.tmp > %s.txt' %(out_file, out_file))
	os.system('rm %s.tmp header.txt' %(out_file))
