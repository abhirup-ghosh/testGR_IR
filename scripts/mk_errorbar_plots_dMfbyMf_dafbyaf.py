#!/home/abhirup/bin/python

def gf(P):
        return filter.gaussian_filter(P, sigma=2.0)

#import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage.filters as filter
#import plotsettings
import bayesian as ba
import time
import os
import os.path

N_bins = 201
P_dMfbyMf_dafbyaf_joint = np.ones((N_bins, N_bins))
s1_dMfbyMf = np.array([])
s1_dafbyaf = np.array([])
s2_dMfbyMf = np.array([])
s2_dafbyaf = np.array([])
mean_dMfbyMf = np.array([])
mean_dafbyaf = np.array([])
median_dMfbyMf = np.array([])
median_dafbyaf = np.array([])
left1_dMfbyMf = np.array([])
left1_dafbyaf = np.array([])
left2_dMfbyMf = np.array([])
left2_dafbyaf = np.array([])
right1_dMfbyMf = np.array([])
right1_dafbyaf = np.array([])
right2_dMfbyMf = np.array([])
right2_dafbyaf = np.array([])


spin_fit_formula_list = ['nonprecspin_Healy2014']
injection_no_list = list(np.linspace(1,50,50))

in_dir_root = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-11-18/uniform_mtot_q_prior_spins/imrtestgr_nonprecspin_Healy2014'

for injection_no in injection_no_list:
        in_dir = '%s/injection_%d/data'%(in_dir_root, injection_no)
        if os.path.isdir(in_dir) == True:

		P_dMfbyMf_dafbyaf = np.loadtxt('%s/P_dMfbyMf_dafbyaf.dat'%in_dir)
	        dMfbyMf_vec = np.loadtxt('%s/dMfbyMf_vec.dat'%in_dir)
        	dafbyaf_vec = np.loadtxt('%s/dafbyaf_vec.dat'%in_dir)

        	dx = np.mean(np.diff(dMfbyMf_vec))
          	dy = np.mean(np.diff(dafbyaf_vec))

          	# Normalization

          	P_dMfbyMf_dafbyaf /= np.sum(P_dMfbyMf_dafbyaf) * dx * dy

          	# Marginalization to one-dimensional joint_posteriors

          	P_dMfbyMf = np.sum(P_dMfbyMf_dafbyaf, axis=0) * dy
          	P_dafbyaf = np.sum(P_dMfbyMf_dafbyaf, axis=1) * dx

		# Calculation of confidence levels

                s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf, 0.68)
                s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf, 0.95)

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

                left1_v2 = min(dafbyaf_vec[np.where(P_dafbyaf>=s1_v2)[0]])
                right1_v2 = max(dafbyaf_vec[np.where(P_dafbyaf>=s1_v2)[0]])

                left2_v2 = min(dafbyaf_vec[np.where(P_dafbyaf>=s2_v2)[0]])
                right2_v2 = max(dafbyaf_vec[np.where(P_dafbyaf>=s2_v2)[0]])

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

	
plt.figure(figsize=(20,8))
plt.subplot(211)
plt.plot(1+np.arange(len(mean_dMfbyMf)), (mean_dMfbyMf), '.', marker='o')
plt.errorbar(1+np.arange(len(mean_dMfbyMf)), np.zeros(len(mean_dMfbyMf)), yerr=(abs(left1_dMfbyMf), abs(right1_dMfbyMf)))
plt.axhline(y=0.)
plt.subplot(212)
plt.plot(1+np.arange(len(mean_dafbyaf)), (mean_dafbyaf), '.', marker='o')
plt.errorbar(1+np.arange(len(mean_dafbyaf)), np.zeros(len(mean_dafbyaf)), yerr=(abs(left1_dafbyaf), abs(right1_dafbyaf)))
plt.axhline(y=0.)
plt.savefig('%s/confidence_plots_casewise.png'%in_dir_root)
