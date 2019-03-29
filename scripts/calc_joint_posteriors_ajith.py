"""

Compute the joint posterior on Delta Mf/Mf and Delta af/af from multiple simulated BBH events

A. Ghosh, P. Ajith, 2015-11-27

 """
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage.filters as filter
import tgrplotsettings
import bayesian as ba
import time
import os
import os.path
import random
from random import shuffle


def gf(P):
	return filter.gaussian_filter(P, sigma=2.0)

def optimal_snr_module(filename):
	data = np.genfromtxt(filename, dtype=None, names=True, usecols=(0,1,2,3,4,5,6))
	var_names = [d[0] for d in data]
	stat_names = data.dtype.names
	optimal_snr = data[var_names.index('optimal_snr')][stat_names.index('mean')+1]	
	return optimal_snr

def calc_conf_intervals_in_1d(P, x):

		# find the value of P corresponding to 68% and 95% confidence heights 
		P_s1 = ba.nsigma_value(P, 0.68)
		P_s2 = ba.nsigma_value(P, 0.95)

		# calculation of condifence edges (values of x corresponding to the height s1 on the two sides) 
		x_s1_l = min(x[np.where(P >= P_s1)[0]])
		x_s1_r = max(x[np.where(P >= P_s1)[0]])

		# calculation of condifence edges (values of x corresponding to the height s2 on the two sides) 
		x_s2_l = min(x[np.where(P >= P_s2)[0]])
		x_s2_r = max(x[np.where(P >= P_s2)[0]])

		return P_s1, P_s2, x_s1_l, x_s1_r, x_s2_l, x_s2_r 

##### MAIN CODE ########
mk_debug_plots = False 
spin_fit_formula = 'nonprecspin_Healy2014'
injection_no_list = np.arange(1, 100)
#random.shuffle(injection_no_list)
SNR_THRESH = 0.
	
in_dir_root = '/home/abhirup/public_html/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-11-18/uniform_mtot_q_prior_spins/imrtestgr_%s_range_-2.0_2.0' %(spin_fit_formula)
out_dir     = '/home/ajith/public_html/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-11-18/uniform_mtot_q_prior_spins/imrtestgr_%s_range_-2.0_2.0_snrthresh%d' %(spin_fit_formula, SNR_THRESH)
mark_size = 2

# create output directory 
os.system('mkdir -p %s'%out_dir) 

# initialize variables 
N_bins = 401
P_dMfbyMf_dafbyaf_joint = np.ones((N_bins, N_bins))
mean_dMfbyMf_joint = np.array([])
mean_dafbyaf_joint = np.array([])
left1_dMfbyMf_joint = np.array([])
left1_dafbyaf_joint = np.array([])
left2_dMfbyMf_joint = np.array([])
left2_dafbyaf_joint = np.array([])
right1_dMfbyMf_joint = np.array([])
right1_dafbyaf_joint = np.array([])
right2_dMfbyMf_joint = np.array([])
right2_dafbyaf_joint = np.array([])

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
insp_snr_arr = np.array([])
ring_snr_arr = np.array([])
imr_snr_arr = np.array([])

plt.figure(figsize=(7,3.8))
plt.subplot2grid((3,8), (1,5), rowspan=2, colspan=3)

for injection_no in injection_no_list:
	in_dir = '%s/injection_%d' %(in_dir_root, injection_no)

	if os.path.isfile('%s/data/P_dMfbyMf_dafbyaf.dat'%in_dir) == True:	

		print '... processing injection ', injection_no
  
		# read the text files containing the optimal snr 
		insp_summary_stats_loc = '%s/lalinf_insp/summary_statistics.dat'%in_dir
		ring_summary_stats_loc = '%s/lalinf_ring/summary_statistics.dat'%in_dir
		imr_summary_stats_loc = '%s/lalinf_imr/summary_statistics.dat'%in_dir

		# read the optimal from each of the runs 
		insp_snr = optimal_snr_module(insp_summary_stats_loc)
		ring_snr = optimal_snr_module(ring_summary_stats_loc)
		imr_snr = optimal_snr_module(imr_summary_stats_loc)

		# apply some threshold 
		if insp_snr > SNR_THRESH and ring_snr > SNR_THRESH and imr_snr > SNR_THRESH:

			# append to the SNR vector (SNRs for each event) 
			insp_snr_arr = np.append(insp_snr_arr, insp_snr)
			ring_snr_arr = np.append(ring_snr_arr, ring_snr)
			imr_snr_arr = np.append(imr_snr_arr, imr_snr)

			# read the posterior data for this event 
			P_dMfbyMf_dafbyaf = np.loadtxt('%s/data/P_dMfbyMf_dafbyaf.dat'%in_dir)
			dMfbyMf_vec = np.loadtxt('%s/data/dMfbyMf_vec.dat'%in_dir)
			dafbyaf_vec = np.loadtxt('%s/data/dafbyaf_vec.dat'%in_dir)

			# Joint Probability distribution computation
			P_dMfbyMf_dafbyaf_joint = P_dMfbyMf_dafbyaf_joint*P_dMfbyMf_dafbyaf
			dx = np.mean(np.diff(dMfbyMf_vec))
			dy = np.mean(np.diff(dafbyaf_vec))
			P_dMfbyMf_dafbyaf_joint /= np.sum(P_dMfbyMf_dafbyaf_joint) * dx * dy	# normalization 

			# Marginalization to one-dimensional joint_posteriors
			P_dMfbyMf_joint = np.sum(P_dMfbyMf_dafbyaf_joint, axis=0) * dy
			P_dafbyaf_joint = np.sum(P_dMfbyMf_dafbyaf_joint, axis=1) * dx
			P_dMfbyMf = np.sum(P_dMfbyMf_dafbyaf, axis=0) * dy
			P_dafbyaf = np.sum(P_dMfbyMf_dafbyaf, axis=1) * dx

			# Calculation of confidence levels in the 2D posterior 
			s1_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dafbyaf_joint, 0.68)
			s2_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dafbyaf_joint, 0.95)
			s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf, 0.68)
			s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf, 0.95)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (joint) 
			s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_conf_intervals_in_1d(P_dMfbyMf_joint, dMfbyMf_vec)
			s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_conf_intervals_in_1d(P_dafbyaf_joint, dafbyaf_vec)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (current injection) 
			s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_conf_intervals_in_1d(P_dMfbyMf, dMfbyMf_vec)
			s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_conf_intervals_in_1d(P_dafbyaf, dafbyaf_vec)

			# append the confidence intervals for plotting (joint posteriors) 
			left1_dMfbyMf_joint = np.append(left1_dMfbyMf_joint, left1_v1_joint)
			right1_dMfbyMf_joint = np.append(right1_dMfbyMf_joint, right1_v1_joint)
			left2_dMfbyMf_joint = np.append(left2_dMfbyMf_joint, left2_v1_joint)
			right2_dMfbyMf_joint = np.append(right2_dMfbyMf_joint, right2_v1_joint)
			left1_dafbyaf_joint = np.append(left1_dafbyaf_joint, left1_v2_joint)
			right1_dafbyaf_joint = np.append(right1_dafbyaf_joint, right1_v2_joint)
			left2_dafbyaf_joint = np.append(left2_dafbyaf_joint, left2_v2_joint)
			right2_dafbyaf_joint = np.append(right2_dafbyaf_joint, right2_v2_joint)

			# compute the mean value of the conf interval (joint posteriors) 
			mean_dMfbyMf_joint = np.append(mean_dMfbyMf_joint, np.average(dMfbyMf_vec, weights=P_dMfbyMf_joint))
			mean_dafbyaf_joint = np.append(mean_dafbyaf_joint, np.average(dafbyaf_vec, weights=P_dafbyaf_joint))

			# append the confidence intervals for plotting (current injection) 
			left1_dMfbyMf = np.append(left1_dMfbyMf, left1_v1)
			right1_dMfbyMf = np.append(right1_dMfbyMf, right1_v1)
			left2_dMfbyMf = np.append(left2_dMfbyMf, left2_v1)
			right2_dMfbyMf = np.append(right2_dMfbyMf, right2_v1)
			left1_dafbyaf = np.append(left1_dafbyaf, left1_v2)
			right1_dafbyaf = np.append(right1_dafbyaf, right1_v2)
			left2_dafbyaf = np.append(left2_dafbyaf, left2_v2)
			right2_dafbyaf = np.append(right2_dafbyaf, right2_v2)

			# compute the mean value of the conf interval (current injection) 
			mean_dMfbyMf = np.append(mean_dMfbyMf, np.average(dMfbyMf_vec, weights=P_dMfbyMf))
			mean_dafbyaf = np.append(mean_dafbyaf, np.average(dafbyaf_vec, weights=P_dafbyaf))

			# make a debug plot 
			if mk_debug_plots == True: 

				plt.figure(figsize=(5,5))
				plt.subplot2grid((3,3), (0,0), colspan=2)
				plt.plot(dMfbyMf_vec, P_dMfbyMf, color='k', lw=1)
				plt.axvline(x=left1_v1, color='c', lw=0.5, ls='-.')
				plt.axvline(x=right1_v1, color='c', lw=0.5, ls='-.')
				plt.axvline(x=left2_v1, color='b', lw=0.5, ls='-.')
				plt.axvline(x=right2_v1, color='b', lw=0.5, ls='-.')
				plt.ylabel('$P(\Delta M_f/M_f)$')
				plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
				plt.pcolormesh(dMfbyMf_vec,dafbyaf_vec,P_dMfbyMf_dafbyaf, cmap='YlOrBr')
				plt.contour(dMfbyMf_vec,dafbyaf_vec, gf(P_dMfbyMf_dafbyaf), levels=(s1_v1v2,s2_v1v2), linewidths=(1,1.5))
				plt.plot(0, 0, 'k+', ms=12, mew=2)
				plt.xlabel('$\Delta M_f/M_f$')
				plt.ylabel('$\Delta a_f/a_f$')
				plt.xlim([-1.,1.])
				plt.ylim([-1.,1.])
				plt.grid()
				plt.subplot2grid((3,3), (1,2), rowspan=2)
				plt.plot(P_dafbyaf, dafbyaf_vec,'k', lw=1)
				plt.axhline(y=left1_v2, color='c', lw=0.5, ls='-.')
				plt.axhline(y=right1_v2, color='c', lw=0.5, ls='-.')
				plt.axhline(y=left2_v2, color='b', lw=0.5, ls='-.')
				plt.axhline(y=right2_v2, color='b', lw=0.5, ls='-.')
				plt.xlabel('$P(\Delta a_f/a_f)$')
				plt.savefig('%s/dMfbyMfdafbyaf_inj%d.png' %(out_dir, injection_no), dpi=300)
				plt.close()

			else: 
				# plot the 2D 1-sigma contour from the current posteriors 
				plt.contour(dMfbyMf_vec,dafbyaf_vec, gf(P_dMfbyMf_dafbyaf), levels=(s1_v1v2,), linewidths=(0.1,), alpha=0.5, colors='#ff7c4c')
				plt.contour(dMfbyMf_vec,dafbyaf_vec, gf(P_dMfbyMf_dafbyaf_joint), levels=(s1_v1v2,), linewidths=(0.1,), alpha=float(injection_no)/len(injection_no_list), colors='r')
				plt.hold(True)

#########################################################################
######## plot the confidence intervals against the event number  ########
#########################################################################

upper_error_dMfbyMf = mean_dMfbyMf - left1_dMfbyMf
lower_error_dMfbyMf = right1_dMfbyMf - mean_dMfbyMf
upper_error_dafbyaf = mean_dafbyaf - left1_dafbyaf
lower_error_dafbyaf = right1_dafbyaf - mean_dafbyaf  


plt.contour(dMfbyMf_vec,dafbyaf_vec,gf(P_dMfbyMf_dafbyaf_joint), levels=(s1_v1v2_joint,s2_v1v2_joint), linewidths=(.3,1), colors='#cc0000')
plt.plot(0, 0, 'k+', mew=0.5)
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta a_f/a_f$')
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta a_f/a_f$')

ax = plt.subplot2grid((3,8), (0,0), colspan=4)
ax.set_xticklabels([])
plt.plot(1+np.arange(len(insp_snr_arr)), insp_snr_arr, 'co-', label='$\\rho_\mathrm{insp}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(ring_snr_arr)), ring_snr_arr, 'ro-', label='$\\rho_\mathrm{post-insp}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(imr_snr_arr)), imr_snr_arr, 'ko-', label='$\\rho_\mathrm{IMR}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(imr_snr_arr)), np.sqrt(insp_snr_arr**2.+ring_snr_arr**2.), label='$(\\rho_\mathrm{insp}^2+\\rho_\mathrm{post-insp}^2)^{1/2}$', lw=1, mew=0)
plt.legend(loc='lower right', fontsize=5, frameon=False)
plt.xlim(1, len(mean_dMfbyMf))
plt.ylim(0,30)
plt.ylabel('mean $\\rho_\mathrm{opt}$')
plt.grid()
ax = plt.subplot2grid((3,8), (1,0), colspan=4)
ax.set_xticklabels([])
plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left1_dMfbyMf_joint, right1_dMfbyMf_joint, color='#ff7c4c', alpha=0.8) #, label='68\%')
plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left2_dMfbyMf_joint, right2_dMfbyMf_joint, color='#ff7c4c', alpha=0.2) #, label='95\%')
plt.plot(1+np.arange(len(mean_dMfbyMf)), (mean_dMfbyMf), '.', marker='o', color='#cc0000', mew=0, ms=mark_size)
plt.errorbar(1+np.arange(len(mean_dMfbyMf)), mean_dMfbyMf, yerr=(upper_error_dMfbyMf, lower_error_dMfbyMf), linestyle='None', alpha=1, color='#cc0000', lw=0.1, capsize=0)
plt.axhline(y=0., color='k', ls='--', lw=0.5)
plt.ylim(-0.4,0.4)
plt.xlim(1, len(mean_dMfbyMf))
plt.ylabel('$\Delta M_f/M_f$')
plt.subplot2grid((3,8), (2,0), colspan=4)
plt.fill_between(1+np.arange(len(mean_dafbyaf_joint)), left1_dafbyaf_joint, right1_dafbyaf_joint, color='#ff7c4c', alpha=0.8) #, label='68\%')
plt.fill_between(1+np.arange(len(mean_dafbyaf_joint)), left2_dafbyaf_joint, right2_dafbyaf_joint, color='#ff7c4c', alpha=0.2) #, label='95\%')
plt.plot(1+np.arange(len(mean_dafbyaf)), (mean_dafbyaf), '.', marker='o', color='#cc0000', mew=0, ms=mark_size)
plt.errorbar(1+np.arange(len(mean_dafbyaf)), mean_dafbyaf, yerr=(upper_error_dafbyaf, lower_error_dafbyaf), linestyle='None', alpha=1, color='#cc0000', lw=0.1, capsize=0)
plt.axhline(y=0., color='k', ls='--', lw=0.5)
plt.ylim(-0.4,0.4)
plt.xlim(1, len(mean_dMfbyMf))
plt.xlabel('Number of events')
plt.ylabel('$\Delta a_f/a_f$')


plt.savefig('%s/confidence_levels_joint_posterior.png' %out_dir, dpi=300)
plt.close()

