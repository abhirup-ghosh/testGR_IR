import matplotlib as mpl
mpl.use('Agg')
import sys, os, numpy as np
import commands
import matplotlib.pyplot as plt
import os
import commands
import glob
import lal

data = np.genfromtxt("./SXS_campaign.dat", names=True, dtype=None)

for idx in range(len(data)):
        print idx, data[idx]['tag']

        tag = data[idx]['tag']

        post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/systematics_error_characterisation/sxs_injection_o1o2_noise/%s'%(tag)
	out_dir = post_loc + '/imrtgr_results'

	m1, m2, s1x, s1y, s1z, s2x, s2y, s2z = data[idx]['m1'],data[idx]['m2'], data[idx]['s1x'], data[idx]['s1y'], data[idx]['s1z'],data[idx]['s2x'], data[idx]['s2y'], data[idx]['s2z']

        s1 = np.sqrt(s1x**2. + s1y**2. + s1z**2.)
        s2 = np.sqrt(s2x**2. + s2y**2. + s2z**2.)
        phi12 = np.arctan(s2y/s2x) - np.arctan(s1y/s1x)
	if np.isnan(phi12) == True:
		phi12 = 0

	# calculate the Kerr ISCO freq 
        Mf, af, f_isco = data[idx]['Mf'],data[idx]['af'], data[idx]['f_isco']


	try:
	    post_samples_i = glob.glob(post_loc + '/inspiral/*/*/*/*/posterior_samples.dat')[0]
	    post_samples_r = glob.glob(post_loc + '/post-inspiral/*/*/*/*/posterior_samples.dat')[0]
	    post_samples_imr = glob.glob(post_loc + '/IMR/*/*/*/*/posterior_samples.dat')[0]

	  

	    if os.path.isfile(post_samples_i) == True and os.path.isfile(post_samples_r) == True and os.path.isfile(post_samples_imr) == True:

		run_cmd = 'python /home/abhirup/opt/lalsuite_master_20190628_e5b55248/libexec/lalinference/imrtgr_imr_consistency_test.py --insp-post=%s --ring-post=%s --imr-post=%s --fit-formula=bbh_average_fits_precessing --out-dir=%s --m1-inj=%f --m2-inj=%f --chi1-inj=%f --chi2-inj=%f --chi1z-inj=%f --chi2z-inj=%f --phi12-inj=%f --insp-fhigh=%f --ring-flow=%f --waveform=SEOBNRv4_ROM --N_bins=201 --dMfbyMf_lim=2 --dchifbychif_lim=2 --use_KDE=0' %(post_samples_i, post_samples_r, post_samples_imr, out_dir, m1, m2, s1, s2, s1z, s2z, phi12, f_isco, f_isco)
		
		print run_cmd 
		os.system(run_cmd)  
		
	except:
	  print 'data not found'
