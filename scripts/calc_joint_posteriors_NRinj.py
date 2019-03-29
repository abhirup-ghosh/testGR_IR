"""

Compute the joint posterior on Delta Mf/Mf and Delta chif/chif from multiple simulated BBH events

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
import commands
import imrtestgr as tgr
import pickle, gzip


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
fit_formula = 'nonprecspin_Healy2014'
injection_no_list = [1]#np.arange(1, 100)
a_vec = np.linspace(0.5, 1, 100)
cat_id = 0 
#random.shuffle(injection_no_list)
SNR_THRESH = 1.
GR_CONF_THRESH = 100.
Q_THRESH = 20.
MTOT_THRESH = 200.
prior_correction = False 
	
in_dir_root = '/home/archis/Work/WaveformSystematics/imrtgr_runs/2015-12-17/SXS-Ossokine'
inj_dir_root = '/home/archis/Work/WaveformSystematics/imrtgr_runs/2015-12-17/SXS-Ossokine'
out_dir = '/home/abhirup/Documents/Work/WaveformSystematics/imrtgr_runs/2015-12-17/SXS-Ossokine/imrtestgr_%s/snrthresh%d_GRthresh_%d_Mthresh_%d_qthresh_%d_catalog%d' %(fit_formula, SNR_THRESH, GR_CONF_THRESH, MTOT_THRESH, Q_THRESH, cat_id)
mark_size = 2

# create output directory 
os.system('mkdir -p %s'%out_dir) 

# initialize variables 
N_bins = 201
P_dMfbyMf_dchifbychif_joint = np.ones((N_bins, N_bins))
mean_dMfbyMf_joint = np.array([])
mean_dchifbychif_joint = np.array([])
left1_dMfbyMf_joint = np.array([])
left1_dchifbychif_joint = np.array([])
left2_dMfbyMf_joint = np.array([])
left2_dchifbychif_joint = np.array([])
right1_dMfbyMf_joint = np.array([])
right1_dchifbychif_joint = np.array([])
right2_dMfbyMf_joint = np.array([])
right2_dchifbychif_joint = np.array([])

mean_dMfbyMf = np.array([])
mean_dchifbychif = np.array([])
left1_dMfbyMf = np.array([])
left1_dchifbychif = np.array([])
left2_dMfbyMf = np.array([])
left2_dchifbychif = np.array([])
right1_dMfbyMf = np.array([])
right1_dchifbychif = np.array([])
right2_dMfbyMf = np.array([])
right2_dchifbychif = np.array([])
insp_snr_arr = np.array([])
ring_snr_arr = np.array([])
imr_snr_arr = np.array([])
gr_conf_arr = np.array([])
inj_no_arr = np.array([])
Mf_inj_arr = np.array([])
chif_inj_arr = np.array([])
m1_inj_arr = np.array([])
m2_inj_arr = np.array([])
chi1_inj_arr = np.array([])
chi2_inj_arr = np.array([])

m1_inj=40.834475743889904
m2_inj=33.263425054310098
chi1_inj=0.326
chi2_inj=-0.558
M_inj = m1_inj+m2_inj
q_inj = m1_inj/m2_inj
Mf_inj, chif_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, fit_formula)

plt.figure(figsize=(7,3.8))
plt.subplot2grid((3,8), (1,5), rowspan=2, colspan=3)

for injection_no in injection_no_list:
  in_dir_list = commands.getoutput('ls -d %s/*/*/*/imrtgr'%in_dir_root).split()
  #np.random.shuffle(in_dir_list)
  for in_dir in in_dir_list:
    if os.path.isfile('%s/data/P_dMfbyMf_dchifbychif.dat.gz'%in_dir) == True:
        print in_dir
        print '... processing injection ', injection_no

	if os.path.isfile('%s/data/P_dMfbyMf_dchifbychif.dat.gz' %in_dir) == False:	
		print '... data not found'
	else: 
		print '... found data', 

		# read the text files containing the optimal snr 
		insp_summary_stats_loc = '%s/lalinf_insp/summary_statistics.dat'%in_dir
		ring_summary_stats_loc = '%s/lalinf_ring/summary_statistics.dat'%in_dir
		imr_summary_stats_loc = '%s/lalinf_imr/summary_statistics.dat'%in_dir

		insp_snr = optimal_snr_module(insp_summary_stats_loc)
		ring_snr = optimal_snr_module(ring_summary_stats_loc)
		imr_snr = optimal_snr_module(imr_summary_stats_loc)

		# read gr confidence value
		gr_conf_level = np.loadtxt('%s/data/GR_confidence.txt'%in_dir)

		# apply some threshold 
		if insp_snr > SNR_THRESH and ring_snr > SNR_THRESH and imr_snr > SNR_THRESH and gr_conf_level < GR_CONF_THRESH/100.:

			inj_no_arr = np.append(inj_no_arr, injection_no)
			Mf_inj_arr = np.append(Mf_inj_arr, Mf_inj)
			chif_inj_arr = np.append(chif_inj_arr, chif_inj)
			m1_inj_arr = np.append(m1_inj_arr, m1_inj)
			m2_inj_arr = np.append(m2_inj_arr, m2_inj)
			chi1_inj_arr = np.append(chi1_inj_arr, chi1_inj)
			chi2_inj_arr = np.append(chi2_inj_arr, chi2_inj)
			gr_conf_arr = np.append(gr_conf_arr, gr_conf_level)

			# append to the SNR vector (SNRs for each event) 
			insp_snr_arr = np.append(insp_snr_arr, insp_snr)
			ring_snr_arr = np.append(ring_snr_arr, ring_snr)
			imr_snr_arr = np.append(imr_snr_arr, imr_snr)

			# read the posterior data for this event 
			P_dMfbyMf_dchifbychif = np.loadtxt('%s/data/P_dMfbyMf_dchifbychif.dat.gz'%in_dir)
			dMfbyMf_vec = np.loadtxt('%s/data/dMfbyMf_vec.dat.gz'%in_dir)
			dchifbychif_vec = np.loadtxt('%s/data/dchifbychif_vec.dat.gz'%in_dir)

			# Calculation of confidence levels in the 2D posterior 
			s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dchifbychif, 0.68)
			s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dchifbychif, 0.95)

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
				plt.pcolormesh(dMfbyMf_vec,dchifbychif_vec,P_dMfbyMf_dchifbychif, cmap='YlOrBr')
				plt.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif), levels=(s1_v1v2,s2_v1v2), linewidths=(1,1.5))
				plt.plot(0, 0, 'k+', ms=12, mew=2)
				plt.xlabel('$\Delta M_f/M_f$')
				plt.ylabel('$\Delta \chi_f/\chi_f$')
				plt.xlim([-1.,1.])
				plt.ylim([-1.,1.])
				plt.grid()
				plt.subplot2grid((3,3), (1,2), rowspan=2)
				plt.plot(P_dchifbychif, dchifbychif_vec,'k', lw=1)
				plt.axhline(y=left1_v2, color='c', lw=0.5, ls='-.')
				plt.axhline(y=right1_v2, color='c', lw=0.5, ls='-.')
				plt.axhline(y=left2_v2, color='b', lw=0.5, ls='-.')
				plt.axhline(y=right2_v2, color='b', lw=0.5, ls='-.')
				plt.xlabel('$P(\Delta \chi_f/\chi_f)$')
				plt.savefig('%s/dMfbyMfdchifbychif_inj%d.png' %(out_dir, injection_no), dpi=300)
				plt.close()

			else: 
				# plot the 2D 1-sigma contour from the current posteriors 
				plt.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif), levels=(s1_v1v2,), linewidths=(0.1,), alpha=0.5, colors='#ff7c4c')

#
				#plt.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif_joint), levels=(s1_v1v2,), linewidths=(0.1,), alpha=float(injection_no)/len(injection_no_list), colors='c')
				#plt.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif_joint), levels=(s1_v1v2_joint,), linewidths=(0.1,), alpha=a_vec[injection_no], colors='c')
				plt.hold(True)

			print 'prior_correction = ', prior_correction

			# for the second event onwards, take out the prior -- the prior is included only once -- for the first event 
			if prior_correction == True: 
				prior_file = '../data/Prior_dMfbyMf_dafbyaf_Mfmin_0.0_Mfmax_1000.0_afmin0.0_afmax1.0'
				f = gzip.open(prior_file+".pklz",'rb')
				P_dMfbyMf_dchifbychif_pr_interp_obj = pickle.load(f)
				P_dMfbyMf_dchifbychif_pr = P_dMfbyMf_dchifbychif_pr_interp_obj(dMfbyMf_vec, dchifbychif_vec)	
				P_dMfbyMf_dchifbychif /= P_dMfbyMf_dchifbychif_pr

			# for the second event onwards, take out the prior -- the prior is included only once -- for the first event 
			prior_correction = True 

 			# removing nans and infinities 
 			P_dMfbyMf_dchifbychif[np.isnan(P_dMfbyMf_dchifbychif)] = 0.
 			P_dMfbyMf_dchifbychif[np.isinf(P_dMfbyMf_dchifbychif)] = 0.

			# replace all zeros in the posterior vector by a small number. 
			# This is to avoid multiplication by zeros while combining the posteiror 
			zidx = np.where(P_dMfbyMf_dchifbychif == 0.) 
			P_dMfbyMf_dchifbychif[zidx] = 1e-6 

			# Joint Probability distribution computation
			P_dMfbyMf_dchifbychif_joint = P_dMfbyMf_dchifbychif_joint*P_dMfbyMf_dchifbychif
			dx = np.mean(np.diff(dMfbyMf_vec))
			dy = np.mean(np.diff(dchifbychif_vec))
			P_dMfbyMf_dchifbychif_joint /= np.sum(P_dMfbyMf_dchifbychif_joint) * dx * dy	# normalization 

			# Marginalization to one-dimensional joint_posteriors
			P_dMfbyMf_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=0) * dy
			P_dchifbychif_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=1) * dx
			P_dMfbyMf = np.sum(P_dMfbyMf_dchifbychif, axis=0) * dy
			P_dchifbychif = np.sum(P_dMfbyMf_dchifbychif, axis=1) * dx

			# Calculation of confidence levels in the 2D posterior 
			s1_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint, 0.68)
			s2_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint, 0.95)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (joint) 
			s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_conf_intervals_in_1d(P_dMfbyMf_joint, dMfbyMf_vec)
			s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_conf_intervals_in_1d(P_dchifbychif_joint, dchifbychif_vec)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (current injection) 
			s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_conf_intervals_in_1d(P_dMfbyMf, dMfbyMf_vec)
			s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_conf_intervals_in_1d(P_dchifbychif, dchifbychif_vec)

			# append the confidence intervals for plotting (joint posteriors) 
			left1_dMfbyMf_joint = np.append(left1_dMfbyMf_joint, left1_v1_joint)
			right1_dMfbyMf_joint = np.append(right1_dMfbyMf_joint, right1_v1_joint)
			left2_dMfbyMf_joint = np.append(left2_dMfbyMf_joint, left2_v1_joint)
			right2_dMfbyMf_joint = np.append(right2_dMfbyMf_joint, right2_v1_joint)
			left1_dchifbychif_joint = np.append(left1_dchifbychif_joint, left1_v2_joint)
			right1_dchifbychif_joint = np.append(right1_dchifbychif_joint, right1_v2_joint)
			left2_dchifbychif_joint = np.append(left2_dchifbychif_joint, left2_v2_joint)
			right2_dchifbychif_joint = np.append(right2_dchifbychif_joint, right2_v2_joint)

			# compute the mean value of the conf interval (joint posteriors) 
			mean_dMfbyMf_joint = np.append(mean_dMfbyMf_joint, np.average(dMfbyMf_vec, weights=P_dMfbyMf_joint))
			mean_dchifbychif_joint = np.append(mean_dchifbychif_joint, np.average(dchifbychif_vec, weights=P_dchifbychif_joint))

			# append the confidence intervals for plotting (current injection) 
			left1_dMfbyMf = np.append(left1_dMfbyMf, left1_v1)
			right1_dMfbyMf = np.append(right1_dMfbyMf, right1_v1)
			left2_dMfbyMf = np.append(left2_dMfbyMf, left2_v1)
			right2_dMfbyMf = np.append(right2_dMfbyMf, right2_v1)
			left1_dchifbychif = np.append(left1_dchifbychif, left1_v2)
			right1_dchifbychif = np.append(right1_dchifbychif, right1_v2)
			left2_dchifbychif = np.append(left2_dchifbychif, left2_v2)
			right2_dchifbychif = np.append(right2_dchifbychif, right2_v2)

			# compute the mean value of the conf interval (current injection) 
			mean_dMfbyMf = np.append(mean_dMfbyMf, np.average(dMfbyMf_vec, weights=P_dMfbyMf))
			mean_dchifbychif = np.append(mean_dchifbychif, np.average(dchifbychif_vec, weights=P_dchifbychif))


#########################################################################
######## plot the confidence intervals against the event number  ########
#########################################################################

upper_error_dMfbyMf = mean_dMfbyMf - left1_dMfbyMf
lower_error_dMfbyMf = right1_dMfbyMf - mean_dMfbyMf
upper_error_dchifbychif = mean_dchifbychif - left1_dchifbychif
lower_error_dchifbychif = right1_dchifbychif - mean_dchifbychif  


plt.contour(dMfbyMf_vec,dchifbychif_vec,gf(P_dMfbyMf_dchifbychif_joint), levels=(s1_v1v2_joint,s2_v1v2_joint), linewidths=(0.5,1.5), colors='#cc0000')
plt.plot(0, 0, 'k+', mew=0.5)
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta \chi_f/\chi_f$')
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta \chi_f/\chi_f$')

ax = plt.subplot2grid((3,8), (0,0), colspan=4)
ax.set_xticklabels([])
plt.plot(1+np.arange(len(gr_conf_arr)), gr_conf_arr*100, 'go-', label='GR confidence level', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(insp_snr_arr)), insp_snr_arr, 'co-', label='$\\rho_\mathrm{insp}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(ring_snr_arr)), ring_snr_arr, 'ro-', label='$\\rho_\mathrm{post-insp}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(imr_snr_arr)), imr_snr_arr, 'ko-', label='$\\rho_\mathrm{IMR}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(imr_snr_arr)), np.sqrt(insp_snr_arr**2.+ring_snr_arr**2.), label='$(\\rho_\mathrm{insp}^2+\\rho_\mathrm{post-insp}^2)^{1/2}$', lw=1, mew=0)
plt.legend(loc='lower right', fontsize=5, frameon=False)
plt.xlim(1, len(mean_dMfbyMf))
plt.ylim(0,100)
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
plt.fill_between(1+np.arange(len(mean_dchifbychif_joint)), left1_dchifbychif_joint, right1_dchifbychif_joint, color='#ff7c4c', alpha=0.8) #, label='68\%')
plt.fill_between(1+np.arange(len(mean_dchifbychif_joint)), left2_dchifbychif_joint, right2_dchifbychif_joint, color='#ff7c4c', alpha=0.2) #, label='95\%')
plt.plot(1+np.arange(len(mean_dchifbychif)), (mean_dchifbychif), '.', marker='o', color='#cc0000', mew=0, ms=mark_size)
plt.errorbar(1+np.arange(len(mean_dchifbychif)), mean_dchifbychif, yerr=(upper_error_dchifbychif, lower_error_dchifbychif), linestyle='None', alpha=1, color='#cc0000', lw=0.1, capsize=0)
plt.axhline(y=0., color='k', ls='--', lw=0.5)
plt.ylim(-0.4,0.4)
plt.xlim(1, len(mean_dMfbyMf))
plt.xlabel('Number of events')
plt.ylabel('$\Delta \chi_f/\chi_f$')
plt.savefig('%s/dMfbyMf_dchifbychif_jointposterior_vs_events.png'%out_dir, dpi=300)
plt.savefig('%s/dMfbyMf_dchifbychif_jointposterior_vs_events.pdf'%out_dir)
plt.close()
plt.savefig('%s/confidence_levels.png'%out_dir, dpi=300)
plt.close()

plt.figure(figsize=(16,10))
plt.subplot(211)
plt.plot(mean_dMfbyMf, gr_conf_arr*100, 'k.', label='GR confidence level', ms=10, lw=1, mew=0)
for (i,(x, y)) in enumerate(zip(mean_dMfbyMf, gr_conf_arr*100)):
    plt.text(x, y, inj_no_arr[i], color="red", fontsize=12)
plt.xlabel('mean dMfbyMf')
plt.ylabel('GR conf')
plt.subplot(212)
plt.plot(mean_dchifbychif, gr_conf_arr*100, 'k.', label='GR confidence level', ms=10, lw=1, mew=0)
for (i,(x, y)) in enumerate(zip(mean_dchifbychif, gr_conf_arr*100)):
    plt.text(x, y, inj_no_arr[i], color="red", fontsize=12)
plt.xlabel('mean dMfbyMf')
plt.ylabel('GR conf')
plt.savefig('%s/mean_dMfbyMf_gr_conf_correlation.png'%out_dir)
plt.close()

q_inj_arr = m1_inj_arr/m2_inj_arr
chi_eff_inj_arr = (m1_inj_arr*chi1_inj_arr + m2_inj_arr*chi2_inj_arr)/(m1_inj_arr+m2_inj_arr)

plt.figure(figsize=(23,10))
plt.subplot(241)
plt.plot(mean_dMfbyMf, gr_conf_arr*100, 'k.', label='GR confidence level', ms=10, lw=1, mew=0)
for (i,(x, y)) in enumerate(zip(mean_dMfbyMf, gr_conf_arr*100)):
	plt.text(x, y, '%d' %inj_no_arr[i], color="red", fontsize=12)
plt.xlabel('mean $\Delta M_f/M_f$')
plt.ylabel('GR conf')
plt.subplot(242)
plt.plot(mean_dchifbychif, gr_conf_arr*100, 'k.', label='GR confidence level', ms=10, lw=1, mew=0)
for (i,(x, y)) in enumerate(zip(mean_dchifbychif, gr_conf_arr*100)):
	plt.text(x, y, '%d' %inj_no_arr[i], color="red", fontsize=12)
plt.xlabel('mean $\Delta \chi_f/\chi_f$')
plt.ylabel('GR conf')

plt.subplot(243)
plt.plot(mean_dMfbyMf, Mf_inj_arr, 'r.', ms=8) 
plt.xlabel('mean $\Delta M_f/M_f$')
plt.ylabel('injected $M_f$')
plt.subplot(244)
plt.plot(mean_dMfbyMf, chif_inj_arr, 'r.', ms=8) 
plt.xlabel('mean $\Delta M_f/M_f$')
plt.ylabel('injected $\chi_f$')

plt.subplot(245)
plt.scatter(chi1_inj_arr, chi2_inj_arr, s=20, c=mean_dMfbyMf, lw=0)
plt.title('color is mean $\Delta M_f/M_f$')
plt.ylabel('injected $\chi_1$')
plt.ylabel('injected $\chi_2$')
plt.colorbar()

plt.subplot(246)
plt.scatter(q_inj_arr, m1_inj_arr+m2_inj_arr, s=20, c=mean_dMfbyMf, lw=0)
plt.title('color is mean $\Delta M_f/M_f$')
plt.xlabel('injected $q$')
plt.ylabel('injected $M$')
plt.colorbar()

plt.subplot(247)
plt.scatter(q_inj_arr, chi1_inj_arr, s=20, c=mean_dMfbyMf, lw=0)
plt.title('color is mean $\Delta M_f/M_f$')
plt.xlabel('injected $q$')
plt.ylabel('injected $\chi_1$')
plt.colorbar()

plt.subplot(248)
plt.scatter(q_inj_arr, chi_eff_inj_arr, s=20, c=mean_dMfbyMf, lw=0)
plt.title('color is mean $\Delta M_f/M_f$')
plt.xlabel('injected $q$')
plt.ylabel('injected $\chi_\mathrm{eff}$')
plt.colorbar()

plt.savefig('%s/mean_dMfbyMf_correlation_plots.png'%out_dir)
plt.close()

# make the P-P plot
px_vec = np.linspace(0, 1, 100)
py_vec = np.zeros(100)
num_inj = float(len(gr_conf_arr))

for ix, px in enumerate(px_vec):
	py_vec[ix] = len(np.where(gr_conf_arr <= px)[0])/num_inj
	
plt.figure(figsize=(5,5))
plt.plot(px_vec, px_vec, 'c', ls='--')
plt.plot(px_vec, py_vec, 'r')
plt.grid()
plt.xlabel('GR confidence threshold')
plt.ylabel('fraction of events with GR confidence below the threshold')
plt.savefig('%s/pp_plots.png' %out_dir, dpi=300)
plt.close()

# plot the 1-sigma and 2-sigma errors against the injections 
inj_vec = 1+np.arange(len(mean_dMfbyMf_joint))
plt.figure(figsize=(9,6))
plt.subplot(221)
plt.semilogy(inj_vec, abs(left1_dMfbyMf_joint-mean_dMfbyMf_joint)*inj_vec**0.5, 'r')
plt.semilogy(inj_vec, abs(left2_dMfbyMf_joint-mean_dMfbyMf_joint)*inj_vec**0.5, 'orange')
plt.ylabel('combined error on $\Delta M_f/M_f$')
plt.xlabel('event number')
plt.subplot(222)
plt.semilogy(inj_vec, abs(left1_dchifbychif_joint-mean_dchifbychif_joint)*inj_vec**0.5, 'r')
plt.semilogy(inj_vec, abs(left2_dchifbychif_joint-mean_dchifbychif_joint)*inj_vec**0.5, 'orange')
plt.xlabel('event number') 
plt.ylabel('combined error on $\Delta \chi_f/\chi_f$')
plt.subplot(223)
plt.plot(inj_vec, mean_dMfbyMf_joint, 'r')
plt.ylabel('combined mean of $\Delta M_f/M_f$')
plt.xlabel('event number')
plt.subplot(224)
plt.plot(inj_vec, mean_dchifbychif_joint, 'r')
plt.xlabel('event number') 
plt.ylabel('combined mean of $\Delta \chi_f/\chi_f$')

plt.savefig('%s/error_vs_numinj.png' %out_dir, dpi=300)


