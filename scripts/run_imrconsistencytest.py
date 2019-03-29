import numpy as np
import os
import glob
import commands

os.system('. ${HOME}/.profile')
home = os.getenv('HOME')
user = commands.getoutput('whoami')
lal_prefix = os.getenv('LAL_PREFIX')
psd_path = '%s/share/lalsimulation'%lal_prefix
glue_location = os.getenv('GLUE_LOCATION')
pylal_location = os.getenv('PYLAL_LOCATION')

#spin_fit_formula, comp_spin_min, comp_spin_max, comp_mass_min, comp_mass_max = 'nospin_Pan2011', 0., 0., 1., 500.
spin_fit_formula, comp_spin_min, comp_spin_max, comp_mass_min, comp_mass_max = 'nonprecspin_Healy2014', -1., 1., 1., 500.
#spin_fit_formula, comp_spin_min, comp_spin_max, comp_mass_min, comp_mass_max = 'nonprecspin_Husa2015', -1., 1., 1., 500.
run_nohup = True

m1_inj, m2_inj, chi1_inj, chi2_inj = 80., 60., 0., 0.

insp_fhigh = 67
ring_flow =  67
root_tag = 'lalsimulation_SEOBNRv2_ROM_DoubleSpinpseudoFourPN_flow30Hz'
waveform = 'SEOBNRv2_ROM_DS'
date = '2016-11-17'

#tag = '%s_%s_%s_fhigh_insp%dHz_flow_ring_%dHz_%s' %(waveform, spin_fit_formula, date, insp_fhigh, ring_flow, root_tag)
tag = '%s_%s_%s_fhigh_insp%dHz_flow_ring_%dHz' %(root_tag, spin_fit_formula, date, insp_fhigh, ring_flow)
#out_dir = '$HOME/public_html/imrtestgr/O1/G211117/%s' %tag
out_dir = '/home/abhirup/Documents/Work/testGR_IR/runs/paper_2_simulations/gaussian_noise/LALAdLIGO/2016-09-25/debugging_injection_m1_80_m2_60/pe/imrtgr/%s_wpc'%tag
os.system('mkdir -p %s'%out_dir)

# location of posterior samples of trigger
post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/paper_2_simulations/gaussian_noise/LALAdLIGO/2016-09-25/debugging_injection_m1_80_m2_60/pe/%s'%root_tag
post_samples_i = glob.glob('%s/inspiral/lalinferencenest/*/*/*/posterior_samples.dat'%post_loc)[0]
post_samples_r = glob.glob('%s/post-inspiral/lalinferencenest/*/*/*/posterior_samples.dat'%post_loc)[0]
post_samples_imr = glob.glob('%s/IMR/lalinferencenest/*/*/*/posterior_samples.dat'%post_loc)[0]

# location of the prior file 
prior_Mfchif_file = '%s/share/lalinference/imrtgr_prior_data/Prior_Mfaf_nonprec_Healy2014_M1-500_aligned.pklz'%lal_prefix
#prior_Mfchif_file = '%s/share/lalinference/imrtgr_prior_data/Prior_Mfaf_nonprec_Healy2014_M1-500_isotropic.pklz'%lal_prefix
prior_corr = True

log_file = '%s/log.txt'%out_dir
err_file = '%s/err.txt'%out_dir

if run_nohup == True:
	if prior_corr == True:
        	run_cmd = 'nohup %s/bin/imrtgr_imr_consistency_test --insp-post=%s --ring-post=%s --imr-post=%s --mf-chif-prior=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401 > %s 2> %s &' %(lal_prefix, post_samples_i, post_samples_r, post_samples_imr, prior_Mfchif_file, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform, log_file, err_file)
	else:
        	run_cmd = 'nohup %s/bin/imrtgr_imr_consistency_test --insp-post=%s --ring-post=%s --imr-post=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401 > %s 2> %s &' %(lal_prefix, post_samples_i, post_samples_r, post_samples_imr, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform, log_file, err_file)
else:
	if prior_corr == True:
        	run_cmd = '%s/bin/imrtgr_imr_consistency_test --insp-post=%s --ring-post=%s --imr-post=%s --mf-chif-prior=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401' %(lal_prefix, post_samples_i, post_samples_r, post_samples_imr, prior_Mfchif_file, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform)
        else:
        	run_cmd = '%s/bin/imrtgr_imr_consistency_test --insp-post=%s --ring-post=%s --imr-post=%s --fit-formula=%s --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=%s --N_bins=401' %(lal_prefix, post_samples_i, post_samples_r, post_samples_imr, spin_fit_formula, out_dir, m1_inj, m2_inj, chi1_inj, chi2_inj, insp_fhigh, ring_flow, waveform)

print run_cmd
os.system(run_cmd)

