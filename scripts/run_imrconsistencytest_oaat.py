import numpy as np
import os
import glob

#spin_fit_formula, comp_spin_min, comp_spin_max, comp_mass_min, comp_mass_max = 'nospin_Pan2011', 0., 0., 1., 500.
spin_fit_formula, comp_spin_min, comp_spin_max, comp_mass_min, comp_mass_max = 'nonprecspin_Healy2014', -1., 1., 1., 500.
#spin_fit_formula, comp_spin_min, comp_spin_max, comp_mass_min, comp_mass_max = 'nonprecspin_Husa2015', -1., 1., 1., 500.
run_nohup = False

test_no = 'G211117'
insp_fhigh = 100
ring_flow = 100
root_tag = 'prod_runs_C01_SpinFix_wo_priorcorr_401bins'
waveform = 'IMRPhenomPv2'#'IMRPhenomB'#'SEOBNRv2_ROM_DS'#'IMRPhenomPv2'
date = '2016-04-20'

tag = '%s_%s_%s_fhigh_insp%dHz_flow_ring_%dHz_%s' %(waveform, spin_fit_formula, date, insp_fhigh, ring_flow, root_tag)
out_dir = '$HOME/public_html/imrtestgr/O1/G211117/%s' %tag
#out_dir = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/IMR_consistency_robustness_tests/imrtestgr/%s'%tag
os.system('mkdir -p %s'%out_dir)

# location of posterior samples of trigger
post_loc = '/home/abhirup/Documents/Work/O1/G211117/20160305_IMRPhenomPv2_C01_nest_SpinFix'
post_samples_i = '%s/inspiral_100Hz/1135136350.65-0/H1L1/posterior_samples.dat'%post_loc
post_samples_r = '%s/ringdown_100Hz/1135136350.65-0/H1L1/posterior_samples.dat'%post_loc
post_samples_imr = '%s/IMR/1135136350.65-0/H1L1/posterior_samples.dat'%post_loc

# location of the prior file 
prior_Mfchif_file = '/home/abhirup/src/lalsuite/lalinference/python/lalinference/imrtgr/Prior_Mfaf_nonprecspin_Healy2014_comp_mass_min1.0_comp_mass_max500.0_comp_spin_min0.0_comp_spin_max1.0_isotropic.pklz'
#prior_Mfchif_file = '/home/abhirup/src/lalsuite/lalinference/python/lalinference/imrtgr/Prior_Mfaf_nonprecspin_Healy2014_comp_mass_min1.0_comp_mass_max500.0_comp_spin_min-1.0_comp_spin_max1.0_aligned.pklz'

log_file = '%s/log.txt'%out_dir
err_file = '%s/err.txt'%out_dir

if run_nohup == True: 
	run_cmd = 'nohup /home/abhirup/opt/lalsuite_imrtgr/bin/imrtgr_imr_consistency_test --insp-post=%s --ring-post=%s --imr-post=%s --mf-chif-prior=%s --fit-formula=%s --out-dir=%s --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins 401 > %s 2> %s &' %(post_samples_i, post_samples_r, post_samples_imr, prior_Mfchif_file, spin_fit_formula, out_dir, insp_fhigh, ring_flow, waveform, log_file, err_file)
else: 
	run_cmd = '/home/abhirup/opt/lalsuite_imrtgr/bin/imrtgr_imr_consistency_test --insp-post=%s --ring-post=%s --imr-post=%s --mf-chif-prior=%s --fit-formula=%s --out-dir=%s --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins 401' %(post_samples_i, post_samples_r, post_samples_imr, prior_Mfchif_file, spin_fit_formula, out_dir, insp_fhigh, ring_flow, waveform)
		
#print run_cmd 
os.system(run_cmd)
