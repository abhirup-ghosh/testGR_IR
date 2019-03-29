#!/usr/bin/env python
"""
script to run the IMR consistency test on simulations as condor jobs

P. Ajith, 2017-03-14
"""

import lal, sys, os, numpy as np
import commands
import matplotlib.pyplot as plt
import os
import imrtestgr as tgr
import nr_fits as nr
import commands
import glob, string

os.system('. ${HOME}/.profile')
lal_prefix = os.getenv('LAL_PREFIX')

# general parameters
run_script = 'python /home/rahul.kashyap/testGR_IR/scripts/imrtgr_imr_consistency_test_final_rahul.py'
prior_corr = False
waveform = 'SEOBNRv4_ROM'
spin_fit_formula = 'nonprecspin_Healy2014'
N_bins = 401

# GR injections 
post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/LIGO_software_injections_H1L1/IMRPPv2_injection_IMRPPv2_recovery/1130030640_1130066791'
out_dir_root = '%s/imrtgr_noMfafprior_401bins'%post_loc

for arg in sys.argv[1:]:
	# select the injection
	inj  = string.atoi(arg, 10)

	print '... processing injection %d' %inj

	sub_dir = 'injection_%04d' %inj

	inj_file = '/home/abhirup/Documents/Work/testGR_IR/runs/LIGO_software_injections_H1L1/IMRPPv2_injection_IMRPPv2_recovery/1130030640_1130066791/IMRPhenomPv2pseudoFourPN_1130030640_1130066791_m_10_80_SNR_50_100_LALAdLIGO_H1L1_20Hz.xml'
	out_dir = '%s/%s' %(out_dir_root, sub_dir)

	try:
		# reading the injection file
		m1_inj_list = map(float, commands.getoutput('ligolw_print %s -c mass1'%inj_file).split('\n'))
		m2_inj_list = map(float, commands.getoutput('ligolw_print %s -c mass2'%inj_file).split('\n'))
		chi1_inj_list = map(float, commands.getoutput('ligolw_print %s -c spin1z'%inj_file).split('\n'))
		chi2_inj_list = map(float, commands.getoutput('ligolw_print %s -c spin2z'%inj_file).split('\n'))
		m1_inj, m2_inj, chi1_inj, chi2_inj = m1_inj_list[inj-1], m2_inj_list[inj-1], chi1_inj_list[inj-1], chi2_inj_list[inj-1]
		
		print '...... read injection file' 

		print '...... read injection file'

		M_inj = m1_inj + m2_inj
		q_inj = m2_inj/m1_inj
		if q_inj < 1.:
			m1_inj, m2_inj = m2_inj, m1_inj
			chi1_inj, chi2_inj = chi2_inj, chi1_inj
			q_inj = 1./q_inj
		# calculate the mass and spin of the final BH
		Mf, af = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, spin_fit_formula)

		print '...... m1 = %2.1f m2 = %2.1f chi1 = %3.2f chi2 = %3.2f Mf = %2.1f af = %3.2f' %(m1_inj, m2_inj, chi1_inj, chi2_inj, Mf, af)

		# calculate the Kerr ISCO freq
		f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)
		print "reached here"
		# calculate the dominant QNM freq
		#f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)
		print "reached after calc_fqnm_dominant"
		insp_fhigh = f_isco_Kerr
		ring_flow = f_isco_Kerr

		post_samples_i = glob.glob('%s/%s/inspiral/*/*/*/*/posterior_samples.dat' %(post_loc, sub_dir))[0]
		post_samples_r = glob.glob('%s/%s/post-inspiral/*/*/*/*/posterior_samples.dat' %(post_loc, sub_dir))[0]
		post_samples_imr = glob.glob('%s/%s/IMR/*/*/*/*/posterior_samples.dat' %(post_loc, sub_dir))[0]

		if os.path.isfile(post_samples_i) == True and os.path.isfile(post_samples_r) == True and os.path.isfile(post_samples_imr) == True:
			print '...... found posterior files'
			if prior_corr == True:
				print '...... running with prior correction'
				run_cmd = '%s --insp-post=%s --ring-post=%s --imr-post=%s --mf-chif-prior=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=%d' %(run_script, post_samples_i, post_samples_r, post_samples_imr, prior_Mfchif_file, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform, N_bins)
			else:
				print '...... running without prior correction'
				run_cmd = '%s --insp-post=%s --ring-post=%s --imr-post=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401' %(run_script, post_samples_i, post_samples_r, post_samples_imr, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform)


		# run command 
		print(run_cmd)
		os.system(run_cmd) 
	except: 
		print '...... injection file not found (%s)' %inj_file
