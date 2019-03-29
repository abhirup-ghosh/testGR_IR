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
import imrtestgr as tgr


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
inj_no_list = np.arange(1, 5000)
#random.shuffle(inj_no_list)
SNR_THRESH = 8.
GR_CONF_THRESH = 100.
Q_THRESH = 20.
MTOT_THRESH = 150.
Nbins = 401
MAX_NUM_EVENTS = 50.
SPIN_THRESH = 1. 
Nbins = 401 
run_tag = 'snrthresh%d_GRthresh_%d_Mthresh_%d_qthresh_%d_chithresh%3.2f' %(SNR_THRESH, GR_CONF_THRESH, MTOT_THRESH, Q_THRESH, SPIN_THRESH)


tag, post_file = 'IMR', 'P_Mfchif_imr.dat.gz'
tag, post_file = 'inspiral', 'P_Mfchif_i.dat.gz'
tag, post_file = 'post-inspiral', 'P_Mfchif_r.dat.gz'
	
inj_type = 'GR'
#inj_type = 'modGR_a2_20'
inj_type = 'IHES_GR'


if inj_type == 'GR':

	# using the more symmetric definition of the new calculation; that avoids the numerical artifacts in the earlier run 
	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final/imrtgr/2017-03-31/altdef/Nbins401/'	

	# using the same defn as above; now usign the same prior as lalinference; ie, uniform in component masses and spins 
	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-03-31/altdef/Nbins401'	
	out_dir = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/GR/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume/imrtgr/2017-03-31/normbymean_gaussiansmoothed401bins_range2_noMfafprior/Nbins%d/%s/' %(Nbins, run_tag)
	
elif inj_type == 'modGR_a2_20':

	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/modGR_a2_20/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-03-31/altdef/Nbins401' # using the same defn as above; now usign the same prior as lalinference; ie, uniform in component masses and spins 
	out_dir  = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/modGR_a2_20/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-30_3det/uniform_compmass_spins_comoving_volume/imrtgr/2017-03-31/normbymean_gaussiansmoothed401bins_range2_noMfafprior/Nbins%d/%s/' %(Nbins, run_tag)

elif inj_type == 'IHES_GR':
	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/GR_20170407/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-03-31/altdef/Nbins401'
	out_dir = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/IHES_GR/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/20170407/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-04-15/Nbins%d/%s/' %(Nbins, run_tag)


# initialize variables 
P_Mfaf_joint = np.ones((Nbins, Nbins))
Mf_bias_vec = np.linspace(-1., 1., Nbins) 
af_bias_vec = np.linspace(-1., 1., Nbins) 
mean_Mf_joint = np.array([])
mean_af_joint = np.array([])
left1_Mf_joint = np.array([])
left1_af_joint = np.array([])
left2_Mf_joint = np.array([])
left2_af_joint = np.array([])
right1_Mf_joint = np.array([])
right1_af_joint = np.array([])
right2_Mf_joint = np.array([])
right2_af_joint = np.array([])

mean_Mf = np.array([])
mean_af = np.array([])
left1_Mf = np.array([])
left1_af = np.array([])
left2_Mf = np.array([])
left2_af = np.array([])
right1_Mf = np.array([])
right1_af = np.array([])
right2_Mf = np.array([])
right2_af = np.array([])
insp_snr_arr = np.array([])
ring_snr_arr = np.array([])
imr_snr_arr = np.array([])

plt.figure(figsize=(10,5))
plt.subplot2grid((3,8), (1,5), rowspan=2, colspan=3)
mark_size = 2

inj_file = '/home/abhirup/Documents/Work/testGR_IR/injections/popsynth_injections_SEOBNRv2_ROM_DoubleSpinthreePointFivePN.txt'
m1_inj_list, m2_inj_list, chi1_inj_list, chi2_inj_list = np.loadtxt(inj_file, usecols=(0,1,6,7), unpack=True)


# create outiput directory 
os.system('mkdir -p %s'%out_dir) 
os.system('cp %s %s' %(__file__, out_dir))
num_successful_events = 0 

for inj_no in inj_no_list:

	# get the params of this injection 
	m1_inj, m2_inj, chi1_inj, chi2_inj = m1_inj_list[inj_no-1], m2_inj_list[inj_no-1], chi1_inj_list[inj_no-1], chi2_inj_list[inj_no-1]

	# ensure that m1 > m2 
	if m1_inj < m2_inj:
		m1_inj, m2_inj = m2_inj, m1_inj
		chi1_inj, chi2_inj = chi2_inj, chi1_inj

	# set subdirectories 
	if inj_type == 'modGR_a2_20' or inj_type == 'IHES_GR':
		in_dir = '%s/IHES_%04d' %(post_loc, inj_no)
		# mod-GR injections are non-spinning 
		chi1_inj, chi2_inj = np.zeros_like(chi1_inj), np.zeros_like(chi2_inj)
	elif inj_type == 'GR':
		in_dir = '%s/injection_%d' %(post_loc, inj_no)
	else:
		raise ValueError('unknown inj_type')

	# calculate the injected values of the final mass/spin
	M_inj = m1_inj+m2_inj
	q_inj = m1_inj/m2_inj
	Mf_inj, chif_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, fit_formula)
	af_inj = chif_inj

	if os.path.isfile('%s/data/%s' %(in_dir, post_file)) != True:	
		#print '... unable to find file', '%s/data/%s' %(in_dir, post_file)
		print '... unable to find file'
	else: 

		print '... processing injection ', inj_no
  
		# read the text files containing the optimal snr 
		insp_summary_stats_loc = '%s/lalinf_insp/summary_statistics.dat'%in_dir
		ring_summary_stats_loc = '%s/lalinf_ring/summary_statistics.dat'%in_dir
		imr_summary_stats_loc = '%s/lalinf_imr/summary_statistics.dat'%in_dir

		# read the optimal from each of the runs 
		insp_snr = optimal_snr_module(insp_summary_stats_loc)
		ring_snr = optimal_snr_module(ring_summary_stats_loc)
		imr_snr = optimal_snr_module(imr_summary_stats_loc)

		# apply some threshold 
		if insp_snr > SNR_THRESH and ring_snr > SNR_THRESH and imr_snr > SNR_THRESH and Mf_inj < 150:

			# append to the SNR vector (SNRs for each event) 
			insp_snr_arr = np.append(insp_snr_arr, insp_snr)
			ring_snr_arr = np.append(ring_snr_arr, ring_snr)
			imr_snr_arr = np.append(imr_snr_arr, imr_snr)

			# read the posterior data for this event 
			P_Mfaf = np.loadtxt('%s/data/%s' %(in_dir, post_file))
			Mf_bins, af_bins = np.loadtxt('%s/data/Mfchif.dat.gz' %in_dir)

			# compute the bias in the estimation of Mf and af for this injection 
			Mf_bias_vec_inj = (Mf_bins[:-1] + Mf_bins[1:])/(2.*Mf_inj) - 1.
			af_bias_vec_inj = (af_bins[:-1] + af_bins[1:])/(2.*af_inj) - 1.

			# interpolate the posterior on to a fixed bias vector 
			P_Mfaf_interp_object = scipy.interpolate.interp2d(Mf_bias_vec_inj, af_bias_vec_inj, P_Mfaf, fill_value=0., bounds_error=False)
			P_Mfaf = P_Mfaf_interp_object(Mf_bias_vec, af_bias_vec)

			# Joint Probability distribution computation
			P_Mfaf_joint = P_Mfaf_joint*P_Mfaf
			dx = np.mean(np.diff(Mf_bias_vec))
			dy = np.mean(np.diff(af_bias_vec))
			P_Mfaf_joint /= np.sum(P_Mfaf_joint) * dx * dy	# normalization 

			# Marginalization to one-dimensional joint_posteriors
			P_Mf_joint = np.sum(P_Mfaf_joint, axis=0) * dy
			P_af_joint = np.sum(P_Mfaf_joint, axis=1) * dx
			P_Mf = np.sum(P_Mfaf, axis=0) * dy
			P_af = np.sum(P_Mfaf, axis=1) * dx

			# Calculation of confidence levels in the 2D posterior 
			s1_v1v2_joint = ba.nsigma_value(P_Mfaf_joint, 0.68)
			s2_v1v2_joint = ba.nsigma_value(P_Mfaf_joint, 0.95)
			s1_v1v2 = ba.nsigma_value(P_Mfaf, 0.68)
			s2_v1v2 = ba.nsigma_value(P_Mfaf, 0.95)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (joint) 
			s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_conf_intervals_in_1d(P_Mf_joint, Mf_bias_vec)
			s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_conf_intervals_in_1d(P_af_joint, af_bias_vec)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (current injection) 
			s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_conf_intervals_in_1d(P_Mf, Mf_bias_vec)
			s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_conf_intervals_in_1d(P_af, af_bias_vec)

			# append the confidence intervals for plotting (joint posteriors) 
			left1_Mf_joint = np.append(left1_Mf_joint, left1_v1_joint)
			right1_Mf_joint = np.append(right1_Mf_joint, right1_v1_joint)
			left2_Mf_joint = np.append(left2_Mf_joint, left2_v1_joint)
			right2_Mf_joint = np.append(right2_Mf_joint, right2_v1_joint)
			left1_af_joint = np.append(left1_af_joint, left1_v2_joint)
			right1_af_joint = np.append(right1_af_joint, right1_v2_joint)
			left2_af_joint = np.append(left2_af_joint, left2_v2_joint)
			right2_af_joint = np.append(right2_af_joint, right2_v2_joint)

			# compute the mean value of the conf interval (joint posteriors) 
			mean_Mf_joint = np.append(mean_Mf_joint, np.average(Mf_bias_vec, weights=P_Mf_joint))
			mean_af_joint = np.append(mean_af_joint, np.average(af_bias_vec, weights=P_af_joint))

			# append the confidence intervals for plotting (current injection) 
			left1_Mf = np.append(left1_Mf, left1_v1)
			right1_Mf = np.append(right1_Mf, right1_v1)
			left2_Mf = np.append(left2_Mf, left2_v1)
			right2_Mf = np.append(right2_Mf, right2_v1)
			left1_af = np.append(left1_af, left1_v2)
			right1_af = np.append(right1_af, right1_v2)
			left2_af = np.append(left2_af, left2_v2)
			right2_af = np.append(right2_af, right2_v2)

			# compute the mean value of the conf interval (current injection) 
			mean_Mf = np.append(mean_Mf, np.average(Mf_bias_vec, weights=P_Mf))
			mean_af = np.append(mean_af, np.average(af_bias_vec, weights=P_af))

			# plot the 2D 1-sigma contour from the current posteriors 
			if num_successful_events < MAX_NUM_EVENTS:
				plt.contour(Mf_bias_vec, af_bias_vec, gf(P_Mfaf), levels=(s1_v1v2,), linewidths=(.5,), alpha=0.5, colors='#ff7c4c')

			# plot the 2D 1-sigma contour from joint posteriors from 10, 25 and 50 events 
			if num_successful_events == 10 or num_successful_events == 25 or num_successful_events == 50:
				plt.contour(Mf_bias_vec, af_bias_vec, gf(P_Mfaf_joint), levels=(s2_v1v2_joint,), linewidths=(0.5,), alpha=num_successful_events/50., colors='#cc0000')
			num_successful_events += 1

#########################################################################
######## plot the confidence intervals against the event number  ########
#########################################################################

upper_error_Mf = mean_Mf - left1_Mf
lower_error_Mf = right1_Mf - mean_Mf
upper_error_af = mean_af - left1_af
lower_error_af = right1_af - mean_af  


plt.contour(Mf_bias_vec,af_bias_vec,gf(P_Mfaf_joint), levels=(s1_v1v2_joint,), linewidths=(1,), colors='k')
plt.plot(0, 0, 'k+', mew=0.5)
plt.xlim(-0.5,0.5)
plt.ylim(-0.5,0.5)
plt.xlabel('$M_f/M_f^\mathrm{inj} - 1$')
plt.ylabel('$\chi_f/\chi_f^\mathrm{inj} - 1$')

ax = plt.subplot2grid((3,8), (0,0), colspan=4)
ax.set_xticklabels([])
plt.plot(1+np.arange(len(insp_snr_arr)), insp_snr_arr, 'co-', label='$\\rho_\mathrm{insp}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(ring_snr_arr)), ring_snr_arr, 'ro-', label='$\\rho_\mathrm{post-insp}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(imr_snr_arr)), imr_snr_arr, 'ko-', label='$\\rho_\mathrm{IMR}$', ms=mark_size, lw=1, mew=0)
plt.plot(1+np.arange(len(imr_snr_arr)), np.sqrt(insp_snr_arr**2.+ring_snr_arr**2.), label='$(\\rho_\mathrm{insp}^2+\\rho_\mathrm{post-insp}^2)^{1/2}$', lw=1, mew=0)
plt.legend(loc='lower right', fontsize=5, frameon=False)
plt.xlim(1, len(mean_Mf))
plt.ylim(0,30)
plt.ylabel('mean $\\rho_\mathrm{opt}$')
plt.grid()
ax = plt.subplot2grid((3,8), (1,0), colspan=4)
ax.set_xticklabels([])
plt.fill_between(1+np.arange(len(mean_Mf_joint)), left1_Mf_joint, right1_Mf_joint, color='#ff7c4c', alpha=0.8) #, label='68\%')
plt.fill_between(1+np.arange(len(mean_Mf_joint)), left2_Mf_joint, right2_Mf_joint, color='#ff7c4c', alpha=0.2) #, label='95\%')
plt.plot(1+np.arange(len(mean_Mf)), (mean_Mf), '.', marker='o', color='#cc0000', mew=0.1, ms=mark_size)
plt.errorbar(1+np.arange(len(mean_Mf)), mean_Mf, yerr=(upper_error_Mf, lower_error_Mf), linestyle='None', alpha=1, color='#cc0000', lw=0.1, capsize=0)
plt.axhline(y=0., color='k', ls='--', lw=0.5)
plt.ylim(-0.3,0.3)
plt.xlim(1, len(mean_Mf))
plt.grid()
plt.ylabel('$M_f/M_f^\mathrm{inj}-1$')
plt.subplot2grid((3,8), (2,0), colspan=4)
plt.fill_between(1+np.arange(len(mean_af_joint)), left1_af_joint, right1_af_joint, color='#ff7c4c', alpha=0.8) #, label='68\%')
plt.fill_between(1+np.arange(len(mean_af_joint)), left2_af_joint, right2_af_joint, color='#ff7c4c', alpha=0.2) #, label='95\%')
plt.plot(1+np.arange(len(mean_af)), (mean_af), '.', marker='o', color='#cc0000', mew=0.1, ms=mark_size)
plt.errorbar(1+np.arange(len(mean_af)), mean_af, yerr=(upper_error_af, lower_error_af), linestyle='None', alpha=1, color='#cc0000', lw=0.1, capsize=0)
plt.axhline(y=0., color='k', ls='--', lw=0.5)
plt.ylim(-0.3,0.3)
plt.xlim(1, len(mean_Mf))
plt.grid()
plt.xlabel('Number of events')
plt.ylabel('$\chi_f/\chi_f^\mathrm{inj}-1$')
plt.tight_layout()


plt.savefig('%s/Mfaf_bias_posteriors_%s.png' %(out_dir, tag), dpi=300)
plt.savefig('%s/Mfaf_bias_posteriors_%s.pdf' %(out_dir, tag))
plt.close()
