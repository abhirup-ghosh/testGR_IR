import sys, os, numpy as np
import commands
import matplotlib.pyplot as plt
import os
import imrtgrutils_final as tgr
import nr_fits as nr
import commands
import glob

os.system('. ${HOME}/.profile')
lal_prefix = os.getenv('LAL_PREFIX')

inj_list = np.linspace(1, 5000, 5000)
spin_fit_formula, comp_spin_min, comp_spin_max, comp_mass_max = 'nonprecspin_Healy2014', -1., 1., 500.
waveform, prior_Mfchif_file = 'SEOBNRv2_ROM_DS', '/home/abhirup/Documents/Work/testGR_IR/data/Prior_Mfaf_nonprecspin_Healy2014_comp_mass_min1.0_comp_mass_max300.00_comp_spin_min-0.98_comp_spin_max1.0_N_sampl1.000000e+07_gaussiansmoothed.pklz'
#waveform, prior_Mfchif_file = 'IMRPhenomPv2', '%s/share/lalinference/imrtgr_prior_data/Prior_Mfaf_nonprec_Healy2014_M1-500_isotropic.pklz'%lal_prefix

#tag = 'GR'
tag = 'modGR'

if tag == 'GR':
    post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume'
elif tag == 'modGR':
    post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/GR'

out_dir_root = '%s/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final'%(post_loc)

inj_file = '/home/abhirup/Documents/Work/testGR_IR/injections/popsynth_injections_SEOBNRv2_ROM_DoubleSpinthreePointFivePN.txt'
m1_inj_list, m2_inj_list, chi1_inj_list, chi2_inj_list = np.loadtxt(inj_file, usecols=(0,1,6,7), unpack=True)

for inj in inj_list:
	inj = int(inj)

	if tag == 'GR':
	   m1_inj, m2_inj, chi1_inj, chi2_inj = m1_inj_list[inj-1], m2_inj_list[inj-1], chi1_inj_list[inj-1], chi2_inj_list[inj-1]
	elif tag == 'modGR':
 	  m1_inj, m2_inj, chi1_inj, chi2_inj = m1_inj_list[inj-1], m2_inj_list[inj-1], 0., 0.

	M_inj = m1_inj + m2_inj
	q_inj = m2_inj/m1_inj

	if q_inj < 1.:
		m1_inj, m2_inj = m2_inj, m1_inj
		chi1_inj, chi2_inj = chi2_inj, chi1_inj
		q_inj = 1./q_inj

	# calculate the mass and spin of the final BH 
        Mf, af = tgr.calc_final_mass_spin(m1_inj, m2_inj, 0., 0., chi1_inj, chi2_inj, 0., 0., 0., spin_fit_formula)
        
	# calculate the Kerr ISCO freq 
        #f_isco_Kerr = nr.calc_isco_freq(af)/(Mf*lal.MTSUN_SI)

        # calculate the dominant QNM freq 
        #f_qnm = nr.calc_fqnm_dominant_mode(af)/(Mf*lal.MTSUN_SI)

        insp_fhigh = 100.#f_isco_Kerr
        ring_flow = 100.#f_isco_Kerr


	print '=== M_inj = %.2f q = %.2f m1 = %.2f m2 = %.2f chi1 = %.2f chi2 = %.2f Mf_inj = %.2f af_inj = %.2f===' %(M_inj, q_inj, m1_inj, m2_inj, chi1_inj, chi2_inj, Mf, af)

	try:
	  if tag == 'GR':
	    post_samples_i = glob.glob('%s/injection_%d/inspiral/*/*/posterior_samples.dat' %(post_loc, inj))[0]
            post_samples_r = glob.glob('%s/injection_%d/post-inspiral/*/*/posterior_samples.dat' %(post_loc, inj))[0]
            post_samples_imr = glob.glob('%s/injection_%d/IMR/*/*/posterior_samples.dat' %(post_loc, inj))[0]

	  elif tag == 'modGR':
	    post_samples_i = glob.glob('%s/IHES_%04d/inspiral/*/*/*/*/posterior_samples.dat' %(post_loc, inj))[0]
            post_samples_r = glob.glob('%s/IHES_%04d/post-inspiral/*/*/*/*/posterior_samples.dat' %(post_loc, inj))[0]
            post_samples_imr = glob.glob('%s/IHES_%04d/IMR/*/*/*/*/posterior_samples.dat' %(post_loc, inj))[0]

	  

	  if os.path.isfile(post_samples_i) == True and os.path.isfile(post_samples_r) == True and os.path.isfile(post_samples_imr) == True:

	    	if tag == 'GR':
		  out_dir = '%s/injection_%d' %(out_dir_root, inj) 
	    	elif tag == 'modGR':
		  out_dir = '%s/IHES_%04d' %(out_dir_root, inj) 
	        os.system('mkdir -p %s' %out_dir)
		log_file = '%s/log.txt' %out_dir
		err_file = '%s/err.txt' %out_dir

		#run_cmd = 'nohup python imrtgr_imr_consistency_test_delta.py --insp-post=%s --ring-post=%s --imr-post=%s --mf-chif-prior=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401 > %s 2> %s &' %(post_samples_i, post_samples_r, post_samples_imr, prior_Mfchif_file, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform, log_file, err_file)
		run_cmd = 'python imrtgr_imr_consistency_test_final.py --insp-post=%s --ring-post=%s --imr-post=%s --mf-chif-prior=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401' %(post_samples_i, post_samples_r, post_samples_imr, prior_Mfchif_file, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform)
		
		print run_cmd 
		os.system(run_cmd)  
		
	except:
	  print 'data not found'
