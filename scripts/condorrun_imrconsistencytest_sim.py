#!/usr/bin/env python
"""
script to run the IMR consistency test on simulations as condor jobs

Abhirup Ghosh, P. Ajith, 2017-03-14
"""

import lal, sys, os, numpy as np
import commands
import matplotlib.pyplot as plt
import os
import commands
import glob, string


signal_data = np.genfromtxt("../notebooks/signal_params_list.dat", names=True, dtype=None)
noise_data = np.genfromtxt("../notebooks/o2_L1H1V1_512_times.dat", names=True, dtype=None)

for arg in sys.argv[1:]:
	# select the injection
	idx = string.atoi(arg, 10)

        post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/systematics_error_characterisation/new_runs/imrppv2_injection_psdmodel_zeronoise/%d'%(noise_data[idx]['o2_L1H1V1_512_geocentric_times'])
        out_dir = post_loc + '/imrtgr_results'

        m1, m2, s1x, s1y, s1z, s2x, s2y, s2z = signal_data[idx]['m1'],signal_data[idx]['m2'], signal_data[idx]['s1x'], signal_data[idx]['s1y'], signal_data[idx]['s1z'],signal_data[idx]['s2x'], signal_data[idx]['s2y'], signal_data[idx]['s2z']

        s1 = np.sqrt(s1x**2. + s1y**2. + s1z**2.)
        s2 = np.sqrt(s2x**2. + s2y**2. + s2z**2.)
        phi12 = np.arctan(s2y/s2x) - np.arctan(s1y/s1x)

        # calculate the Kerr ISCO freq 
        Mf, af, f_isco = signal_data[idx]['Mf'],signal_data[idx]['af'], signal_data[idx]['f_isco_Kerr']


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
