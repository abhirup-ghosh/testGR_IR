def gf(P):
        return filter.gaussian_filter(P, sigma=2.0)

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage.filters as filter
import plotsettings
import bayesian as ba
import time
import os
import os.path

M_inj_list = [50., 75., 100.]
q_inj_list = [1., 4.]
chi1_inj_list = [-0.75, 0., 0.75]
chi2_inj_list = [-0.75, 0., 0.75]

spin_fit_formula_list = ['nospin_Pan2011', 'nonprecspin_Healy2014']
realization_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20]

for spin_fit_formula in spin_fit_formula_list:
  for M_inj in M_inj_list:
    for q_inj in q_inj_list:
      for (chi1_inj, chi2_inj) in zip(chi1_inj_list, chi2_inj_list):
  	
	out_dir = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LIGO-P1200087-v18-aLIGO_EARLY_HIGH/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-10-24/joint_posteriors_%s/%s_%s_%s_%s'%(spin_fit_formula, M_inj, q_inj, chi1_inj, chi2_inj)

	os.system('mkdir -p %s'%(out_dir))

	N_bins = 201
	P_dMfbyMf_dafbyaf_joint = np.ones((N_bins, N_bins))
	#P_dMfbyMf_joint = np.ones((N_bins, N_bins))
	#P_dafbyaf_joint = np.ones((N_bins, N_bins))
	s1_dMfbyMf = np.array([])
	s1_dafbyaf = np.array([])
	s2_dMfbyMf = np.array([])
	s2_dafbyaf = np.array([])
	mean_dMfbyMf = np.array([])
	mean_dafbyaf = np.array([])
	left1_dMfbyMf = np.array([])
	left1_dafbyaf = np.array([])
	left2_dMfbyMf = np.array([])
	left2_dafbyaf = np.array([])
	right1_dMfbyMf = np.array([])
	right1_dafbyaf = np.array([])
	right2_dMfbyMf = np.array([])
	right2_dafbyaf = np.array([])

	plt.figure(figsize=(16,3))
        for realization in realization_list:	

		in_dir = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LIGO-P1200087-v18-aLIGO_EARLY_HIGH/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-10-24/realization_%d/imrtestgr_%s/M_inj%3.2f_q_inj%3.2f_chi1_inj%.2f_chi2_inj%.2f/data' %(realization, spin_fit_formula, M_inj, q_inj, chi1_inj, chi2_inj)
		
		P_dMfbyMf_dafbyaf = np.loadtxt('%s/P_dMfbyMf_dafbyaf.dat'%in_dir)
	        dMfbyMf_vec = np.loadtxt('%s/dMfbyMf_vec.dat'%in_dir)
        	dafbyaf_vec = np.loadtxt('%s/dafbyaf_vec.dat'%in_dir)

        	dx = np.mean(np.diff(dMfbyMf_vec))
        	dy = np.mean(np.diff(dafbyaf_vec))

        	P_dMfbyMf_dafbyaf_joint = P_dMfbyMf_dafbyaf_joint*P_dMfbyMf_dafbyaf

		# Normalization

	        P_dMfbyMf_dafbyaf_joint /= np.sum(P_dMfbyMf_dafbyaf_joint) * dx * dy

		# Marginalization to one-dimensional joint_posteriors

        	P_dMfbyMf = np.sum(P_dMfbyMf_dafbyaf_joint, axis=0) * dy
        	P_dafbyaf = np.sum(P_dMfbyMf_dafbyaf_joint, axis=1) * dx

        	# Calculation of confidence levels

        	s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf_joint, 0.68)
        	s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf_joint, 0.95)

        	s1_v1 = ba.nsigma_value(P_dMfbyMf, 0.68)
        	s2_v1 = ba.nsigma_value(P_dMfbyMf, 0.95)

        	s1_v2 = ba.nsigma_value(P_dafbyaf, 0.68)
        	s2_v2 = ba.nsigma_value(P_dafbyaf, 0.95)

        	# Calculation of confidence regions from levels

        	r1_v1 = len(np.where(P_dMfbyMf>s1_v1)[0]) * dx
        	r2_v1 = len(np.where(P_dMfbyMf>s2_v1)[0]) * dx

        	r1_v2 = len(np.where(P_dafbyaf>s1_v2)[0]) * dy
        	r2_v2 = len(np.where(P_dafbyaf>s2_v2)[0]) * dy

        	s1_dMfbyMf = np.append(s1_dMfbyMf, r1_v1)
        	s2_dMfbyMf = np.append(s2_dMfbyMf, r2_v1)

        	s1_dafbyaf = np.append(s1_dafbyaf, r1_v2)
        	s2_dafbyaf = np.append(s2_dafbyaf, r2_v2)

		# Calculation of condifence edges

        	left1_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])
        	right1_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])

        	left2_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])
        	right2_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])

	        left1_v2 = min(dafbyaf_vec[np.where(P_dafbyaf>s1_v2)[0]])
        	right1_v2 = max(dafbyaf_vec[np.where(P_dafbyaf>s1_v2)[0]])

	        left2_v2 = min(dafbyaf_vec[np.where(P_dafbyaf>s2_v2)[0]])
        	right2_v2 = max(dafbyaf_vec[np.where(P_dafbyaf>s2_v2)[0]])

	        left1_dMfbyMf = np.append(left1_dMfbyMf, left1_v1)
        	right1_dMfbyMf = np.append(right1_dMfbyMf, right1_v1)

	        left2_dMfbyMf = np.append(left2_dMfbyMf, left2_v1)
        	right2_dMfbyMf = np.append(right2_dMfbyMf, right2_v1)

	        left1_dafbyaf = np.append(left1_dafbyaf, left1_v2)
        	right1_dafbyaf = np.append(right1_dafbyaf, right1_v2)

	        left2_dafbyaf = np.append(left2_dafbyaf, left2_v2)
        	right2_dafbyaf = np.append(right2_dafbyaf, right2_v2)

	        mean_dMfbyMf = np.append(mean_dMfbyMf, np.average(dMfbyMf_vec, weights=P_dMfbyMf))
        	mean_dafbyaf = np.append(mean_dafbyaf, np.average(dafbyaf_vec, weights=P_dafbyaf))

	plt.subplot(121)
	plt.plot(1+np.arange(len(mean_dMfbyMf)), mean_dMfbyMf, 'co-')
	plt.fill_between(1+np.arange(len(mean_dMfbyMf)), left1_dMfbyMf, right1_dMfbyMf, color='c', alpha=0.4) #, label='68\%')
	plt.fill_between(1+np.arange(len(mean_dMfbyMf)), left2_dMfbyMf, right2_dMfbyMf, color='c', alpha=0.1) #, label='95\%')
	plt.axhline(y=0., color='k', ls='--', lw=0.5)
	plt.xlabel('Number of events')
	plt.ylabel('$\Delta M_f/M_f$')

	plt.subplot(122)
	plt.plot(1+np.arange(len(mean_dafbyaf)), mean_dafbyaf, 'co-')
	plt.fill_between(1+np.arange(len(mean_dafbyaf)), left1_dafbyaf, right1_dafbyaf, color='c', alpha=0.4) #, label='68\%')
	plt.fill_between(1+np.arange(len(mean_dafbyaf)), left2_dafbyaf, right2_dafbyaf, color='c', alpha=0.1) #, label='95\%')
	plt.axhline(y=0., color='k', ls='--', lw=0.5)
	plt.xlabel('Number of events')
	plt.ylabel('$\Delta a_f/a_f$')
	plt.suptitle('M = %2.1f q = %2.1f chi1 = %.2f chi2 = %.2f' %(M_inj, q_inj, chi1_inj, chi2_inj))
	plt.savefig('%s/confidence_levels.png'%out_dir, dpi=300)
