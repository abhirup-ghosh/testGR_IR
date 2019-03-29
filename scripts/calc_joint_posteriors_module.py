"""

Compute the joint posterior on Delta Mf/Mf and Delta af/af from multiple simulated BBH events

Abhirup Ghosh, P. Ajith, 2015-11-27

 """
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage.filters as filter
import bayesian as ba
import time
import os
import os.path
import random
from random import shuffle
import commands
import imrtestgr as tgr
import pickle, gzip
from scipy import interpolate
from random import shuffle
import tgrplotsettings

# Module for confidence calculations
class confidence(object):
  def __init__(self, counts):
    # Sort in descending order in frequency
    self.counts_sorted = np.sort(counts.flatten())[::-1]
    # Get a normalized cumulative distribution from the mode
    self.norm_cumsum_counts_sorted = np.cumsum(self.counts_sorted) / np.sum(counts)
    # Set interpolations between heights, bins and levels
    self._set_interp()
  def _set_interp(self):
    self._length = len(self.counts_sorted)
    # height from index
    self._height_from_idx = interpolate.interp1d(np.arange(self._length), self.counts_sorted, bounds_error=False, fill_value=0.)
    # index from height
    self._idx_from_height = interpolate.interp1d(self.counts_sorted[::-1], np.arange(self._length)[::-1], bounds_error=False, fill_value=self._length)
    # level from index
    self._level_from_idx = interpolate.interp1d(np.arange(self._length), self.norm_cumsum_counts_sorted, bounds_error=False, fill_value=1.)
    # index from level
    self._idx_from_level = interpolate.interp1d(self.norm_cumsum_counts_sorted, np.arange(self._length), bounds_error=False, fill_value=self._length)
  def level_from_height(self, height):
    return self._level_from_idx(self._idx_from_height(height))
  def height_from_level(self, level):
    return self._height_from_idx(self._idx_from_level(level))

def gf(P):
	return filter.gaussian_filter(P, sigma=1.0)

def injected_medianrecovered_values(filename):
	data = np.genfromtxt(filename, dtype=None, names=True, usecols=(0,1,2,3,4,5,6))
	var_names = [d[0] for d in data]
	stat_names = data.dtype.names
	optimal_snr = data[var_names.index('optimal_snr')][stat_names.index('median')]
	h1_optimal_snr = data[var_names.index('h1_optimal_snr')][stat_names.index('median')]
	l1_optimal_snr = data[var_names.index('l1_optimal_snr')][stat_names.index('median')]
	v1_optimal_snr = 0.#data[var_names.index('v1_optimal_snr')][stat_names.index('median')]
	m1_rec = data[var_names.index('m1')][stat_names.index('median')]
	m2_rec = data[var_names.index('m2')][stat_names.index('median')]
	a1z_rec = data[var_names.index('a1z')][stat_names.index('median')]
	a2z_rec = data[var_names.index('a2z')][stat_names.index('median')]
	Mf_rec = 0.#data[var_names.index('mf')][stat_names.index('median')]
	af_rec = data[var_names.index('af')][stat_names.index('median')]
	m1_inj = data[var_names.index('m1')][stat_names.index('median')+1]
        m2_inj = data[var_names.index('m2')][stat_names.index('median')+1]
        a1z_inj = data[var_names.index('a1z')][stat_names.index('median')+1]
        a2z_inj = data[var_names.index('a2z')][stat_names.index('median')+1]
        Mf_inj = 0.#data[var_names.index('mf')][stat_names.index('median')+1]
        af_inj = data[var_names.index('af')][stat_names.index('median')+1]
	return optimal_snr, h1_optimal_snr, l1_optimal_snr, v1_optimal_snr, m1_rec, m2_rec, a1z_rec, a2z_rec, Mf_rec, af_rec, m1_inj, m2_inj, a1z_inj, a2z_inj, Mf_inj, af_inj

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


def calc_joint_posterior_module(post_loc_list, out_dir, prior_file, snr_thresh, GR_conf_thresh, q_thresh, mtot_thresh, Nbins, MAX_NUM_EVENTS):

  # create epsilon and sigma vectors. The posteriors will be interpolated to these points
  epsilon_vec = np.linspace(-2, 2, Nbins)
  sigma_vec = np.linspace(-2, 2, Nbins)

  # create outiput directory
  os.system('mkdir -p %s'%out_dir)
  os.system('cp %s %s' %(__file__, out_dir))
  os.system('cp %s %s'%(prior_file+".pklz", out_dir))

  # read the prior file and evaluate the prior
  f = gzip.open(prior_file+".pklz",'rb')
  P_dMfbyMf_dafbyaf_pr_interp_obj = pickle.load(f)

  for cat_id in range(3):
 	print '# cat_id = %d' %cat_id
 	num_successful_events = 1

 	# initialize variables
 	P_dMfbyMf_dafbyaf_joint = np.ones((Nbins, Nbins))
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
 	gr_conf_arr = np.array([])
 	gr_conf_eps_arr = np.array([])
 	gr_conf_sigma_arr = np.array([])
 	Mf_inj_arr = np.array([])
 	af_inj_arr = np.array([])
 	m1_inj_arr = np.array([])
	m2_inj_arr = np.array([])
 	a1z_inj_arr = np.array([])
 	a2z_inj_arr = np.array([])

 	fig1 = plt.figure(figsize=(4.6, 4.5))
 	ax1 = fig1.add_subplot(111)
 	left, bottom, width, height = [0.32, 0.7, 0.2, 0.2]
 	ax2 = plt.axes([left, bottom, width, height])

 	fig2 = plt.figure(figsize=(8, 8))
 	f2ax1 = fig2.add_subplot(221)
 	f2ax2 = fig2.add_subplot(222)
 	f2ax3 = fig2.add_subplot(223)
 	f2ax4 = fig2.add_subplot(224)

 	# shuffle the data to crate different catalogs
 	np.random.seed(cat_id)
 	np.random.shuffle(post_loc_list)

        for post_loc in post_loc_list:

            # read the text files containing the optimal snr
                insp_summary_stats_loc = '%s/lalinf_insp/summary_statistics.dat'%post_loc
                ring_summary_stats_loc = '%s/lalinf_ring/summary_statistics.dat'%post_loc
                imr_summary_stats_loc = '%s/lalinf_imr/summary_statistics.dat'%post_loc

                # read the optimal from each of the runs
                insp_snr, insp_h1_snr, insp_l1_snr, insp_v1_snr, insp_m1_rec, insp_m2_rec, insp_a1_rec, insp_a2_rec, insp_Mf_rec, insp_af_rec, insp_m1_inj, insp_m2_inj, insp_a1z_inj, insp_a2z_inj, insp_Mf_inj, insp_af_inj = injected_medianrecovered_values(insp_summary_stats_loc)
                ring_snr, ring_h1_snr, ring_l1_snr, ring_v1_snr, ring_m1_rec, ring_m2_rec, ring_a1_rec, ring_a2_rec, ring_Mf_rec, ring_af_rec, ring_m1_inj, ring_m2_inj, ring_a1z_inj, ring_a2z_inj, ring_Mf_inj, ring_af_inj = injected_medianrecovered_values(ring_summary_stats_loc)
                imr_snr, imr_h1_snr, imr_l1_snr, imr_v1_snr, imr_m1_rec, imr_m2_rec, imr_a1_rec, imr_a2_rec, imr_Mf_rec, imr_af_rec, imr_m1_inj, imr_m2_inj, imr_a1z_inj, imr_a2z_inj, imr_Mf_inj, imr_af_inj = injected_medianrecovered_values(imr_summary_stats_loc)


		m1_inj, m2_inj = imr_m1_inj, imr_m2_inj
		q_inj = m1_inj/m2_inj
		if q_inj < 1.:
			q_inj = 1./q_inj
		M_inj = m1_inj + m2_inj


                # read gr confidence value
                gr_conf_level = np.loadtxt('%s/data/GR_confidence.txt' %post_loc)

		# apply some threshold
                if insp_snr < snr_thresh or ring_snr < snr_thresh or q_inj > q_thresh or M_inj > mtot_thresh or gr_conf_level > GR_conf_thresh/100.:
			print post_loc + ' rejected ...'
		else:

                        Mf_inj_arr = np.append(Mf_inj_arr, imr_Mf_inj)
                        af_inj_arr = np.append(af_inj_arr, imr_af_inj)
                        m1_inj_arr = np.append(m1_inj_arr, imr_m1_inj)
                        m2_inj_arr = np.append(m2_inj_arr, imr_m2_inj)
                        a1z_inj_arr = np.append(a1z_inj_arr, imr_a1z_inj)
                        a2z_inj_arr = np.append(a2z_inj_arr, imr_a2z_inj)

			# read the posterior data for this event
                        P_dMfbyMf_dafbyaf_tmp = np.loadtxt('%s/data/P_dMfbyMf_dchifbychif.dat.gz'%post_loc)
                        epsilon_vec_tmp = np.loadtxt('%s/data/dMfbyMf_vec.dat.gz'%post_loc)
                        sigma_vec_tmp = np.loadtxt('%s/data/dchifbychif_vec.dat.gz'%post_loc)

			# interpolate the posterior to a fixed grid
                        P_dMfbyMf_dafbyaf_obj = scipy.interpolate.interp2d(epsilon_vec_tmp, sigma_vec_tmp, P_dMfbyMf_dafbyaf_tmp, fill_value=0., bounds_error=False)
                        P_dMfbyMf_dafbyaf = P_dMfbyMf_dafbyaf_obj(epsilon_vec, sigma_vec)

			print '...... read posterior files'


			# evaluate the prior at these bins
                        P_dMfbyMf_dafbyaf_pr = P_dMfbyMf_dafbyaf_pr_interp_obj(epsilon_vec, sigma_vec)

                        # Calculation of confidence levels in the 2D posterior
                        s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf, 0.68)
                        s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dafbyaf, 0.95)

                        # plot the 2D 1-sigma contour from the current posteriors
                        if num_successful_events < MAX_NUM_EVENTS:
                                ax1.contour(epsilon_vec,sigma_vec, gf(P_dMfbyMf_dafbyaf), levels=(s1_v1v2,), linewidths=(.5,), alpha=0.5, colors='#ff7c4c')


			# compute the confidence region corresponding to the GR value (delta_Mf/Mf = 0, delta_af/af = 0).
                        # the 'confidence' class is defined on top of this script. This is computed using a flat prior in delta_Mf/Mf, delta_af/af
                        P_dMfbyMf_dafbyaf_flatprior = P_dMfbyMf_dafbyaf/P_dMfbyMf_dafbyaf_pr
                        conf_v1v2 = confidence(P_dMfbyMf_dafbyaf)
                        gr_height = P_dMfbyMf_dafbyaf[np.argmin(abs(epsilon_vec)), np.argmin(abs(sigma_vec))] # taking value closest to (0,0)
                        gr_conf_level = conf_v1v2.level_from_height(gr_height)
                        gr_conf_arr = np.append(gr_conf_arr, gr_conf_level)

                        # removing nans and infinities
                        P_dMfbyMf_dafbyaf[np.isnan(P_dMfbyMf_dafbyaf)] = 0.
                        P_dMfbyMf_dafbyaf[np.isinf(P_dMfbyMf_dafbyaf)] = 0.

                        # replace all zeros in the posterior vector by a small number.
                        # This is to avoid multiplication by zeros while combining the posteiror
                        zidx = np.where(P_dMfbyMf_dafbyaf == 0.)
                        P_dMfbyMf_dafbyaf[zidx] = 1e-5

                        # Joint Probability distribution computation -- from second event onwards, divide the new posterior by the prior so that we are multiplying the
                        # the previous posterior by the current likelihood
                        if num_successful_events == 1:
                                P_dMfbyMf_dafbyaf_joint = P_dMfbyMf_dafbyaf
                        elif num_successful_events > 1:
                                P_dMfbyMf_dafbyaf_joint = P_dMfbyMf_dafbyaf_joint*P_dMfbyMf_dafbyaf/P_dMfbyMf_dafbyaf_pr
                        else:
                                raise ValueError('num_successful_events seems to be less than one!')

                        dx = np.mean(np.diff(epsilon_vec))
                        dy = np.mean(np.diff(sigma_vec))
                        P_dMfbyMf_dafbyaf_joint /= np.sum(P_dMfbyMf_dafbyaf_joint) * dx * dy    # normalization

                        print '...... computed combined posteriors'

			# Marginalization to one-dimensional joint_posteriors
                        P_dMfbyMf_joint = np.sum(P_dMfbyMf_dafbyaf_joint, axis=0) * dy
                        P_dafbyaf_joint = np.sum(P_dMfbyMf_dafbyaf_joint, axis=1) * dx
                        P_dMfbyMf = np.sum(P_dMfbyMf_dafbyaf, axis=0) * dy
                        P_dafbyaf = np.sum(P_dMfbyMf_dafbyaf, axis=1) * dx

                        # calculate hte credible levels in the 1D posterior
                        conf_eps = confidence(P_dMfbyMf)
                        gr_height_eps = P_dMfbyMf[np.argmin(abs(epsilon_vec))]
                        gr_conf_level_eps = conf_eps.level_from_height(gr_height_eps)
                        gr_conf_eps_arr = np.append(gr_conf_eps_arr, gr_conf_level_eps)

                        conf_sigma = confidence(P_dafbyaf)
                        gr_height_sigma = P_dafbyaf[np.argmin(abs(sigma_vec))]
                        gr_conf_level_sigma = conf_sigma.level_from_height(gr_height_sigma)
                        gr_conf_sigma_arr = np.append(gr_conf_sigma_arr, gr_conf_level_sigma)

                        # plot the 1D marginalized posterior on dMfbyMf and dafbyaf
                        f2ax1.plot(epsilon_vec, P_dMfbyMf, linewidth=.5, alpha=0.5, color='#ff7c4c')
                        f2ax2.plot(sigma_vec, P_dafbyaf, linewidth=.5, alpha=0.5, color='#ff7c4c')

                        # Calculation of credible intervals in the 2D posterior
                        s1_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dafbyaf_joint, 0.68)
                        s2_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dafbyaf_joint, 0.95)

                        # calcualte the confidence levels and intervals in the marginalized 1d posteriors (joint)
                        s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_conf_intervals_in_1d(P_dMfbyMf_joint, epsilon_vec)
                        s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_conf_intervals_in_1d(P_dafbyaf_joint, sigma_vec)

                        # calcualte the confidence levels and intervals in the marginalized 1d posteriors (current injection)
                        s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_conf_intervals_in_1d(P_dMfbyMf, epsilon_vec)
                        s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_conf_intervals_in_1d(P_dafbyaf, sigma_vec)


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
                        mean_dMfbyMf_joint = np.append(mean_dMfbyMf_joint, np.average(epsilon_vec, weights=P_dMfbyMf_joint))
                        mean_dafbyaf_joint = np.append(mean_dafbyaf_joint, np.average(sigma_vec, weights=P_dafbyaf_joint))

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
                        mean_dMfbyMf = np.append(mean_dMfbyMf, np.average(epsilon_vec, weights=P_dMfbyMf))
                        mean_dafbyaf = np.append(mean_dafbyaf, np.average(sigma_vec, weights=P_dafbyaf))

			# plot the combined posteriors from 10, 25, 50 events
                        if num_successful_events == int(MAX_NUM_EVENTS/2) or num_successful_events == int(MAX_NUM_EVENTS/5) or num_successful_events == int(MAX_NUM_EVENTS):

                                ax1.contour(epsilon_vec,sigma_vec, P_dMfbyMf_dafbyaf_joint, levels=(s1_v1v2_joint,), linewidths=(1.5,), alpha=num_successful_events/MAX_NUM_EVENTS, colors='#cc0000')
                                ax2.contour(epsilon_vec,sigma_vec, P_dMfbyMf_dafbyaf_joint, levels=(s1_v1v2_joint,), linewidths=(1.5,), alpha=num_successful_events/MAX_NUM_EVENTS, colors='#cc0000')

                                # plot the 1D marginalized posterior on dMfbyMf and dafbyaf
                                f2ax3.plot(epsilon_vec, P_dMfbyMf_joint, linewidth=1.5, alpha=num_successful_events/MAX_NUM_EVENTS, color='#cc0000')
                                f2ax4.plot(sigma_vec, P_dafbyaf_joint, linewidth=1.5, alpha=num_successful_events/MAX_NUM_EVENTS, color='#cc0000')

                        num_successful_events += 1


	#########################################################################
   	######## plot the confidence intervals against the event number  ########
 	#########################################################################

 	upper_error_dMfbyMf = mean_dMfbyMf - left1_dMfbyMf
 	lower_error_dMfbyMf = right1_dMfbyMf - mean_dMfbyMf
 	upper_error_dafbyaf = mean_dafbyaf - left1_dafbyaf
 	lower_error_dafbyaf = right1_dafbyaf - mean_dafbyaf


	# plot 2D posteriors
 	ax1.plot(0, 0, 'k+', mew=1.5, ms=12, alpha=1)
 	ax1.set_xlabel('$\Delta M_f/\\bar{M_f}$', fontsize='xx-large', labelpad=10)
 	ax1.set_ylabel('$\Delta a_f/\\bar{a_f}$', fontsize='xx-large', labelpad=10)
 	ax1.set_xticks([-2, -1, 0, 1, 2])
 	ax1.set_yticks([-2, -1, 0, 1, 2])

 	ax2.plot(0, 0, 'k+', mew=1.5, ms=12, alpha=1)
        ax2.set_xlim([-0.1, 0.1]); ax2.set_ylim([-0.1, 0.1])
 	fig1.tight_layout()
 	fig1.savefig('%s/P_epsilon_sigma_2D_cat%d.pdf' %(out_dir, cat_id))
 	fig1.savefig('%s/P_epsilon_sigma_2D_cat%d.png' %(out_dir, cat_id), dpi=200)
 	plt.close()

 	# plot 1D posteriors
 	f2ax1.set_xlabel('$\Delta M_f/\\bar{M_f}$', fontsize='large', labelpad=10)
 	f2ax1.set_ylabel('$P(\Delta M_f/\\bar{M_f})$', fontsize='large', labelpad=10)
 	f2ax1.set_xlim(min(epsilon_vec), max(epsilon_vec))
 	f2ax3.plot(epsilon_vec, P_dMfbyMf_joint, linewidth=.5, alpha=1, color='k')
 	f2ax3.axvline(x=0)
 	f2ax3.set_xlabel('$\Delta M_f/\\bar{M_f}$', fontsize='large', labelpad=10)
 	f2ax3.set_ylabel('$P(\Delta M_f/\\bar{M_f})$', fontsize='large', labelpad=10)
 	f2ax3.set_xlim(-0.5, 0.5)
 	f2ax2.set_xlabel('$\Delta a_f/\\bar{a_f}$', fontsize='large', labelpad=10)
 	f2ax2.set_ylabel('$P(\Delta a_f/\\bar{a_f})$', fontsize='large', labelpad=10)
 	f2ax2.set_xlim(min(sigma_vec), max(sigma_vec))
 	f2ax4.plot(sigma_vec, P_dafbyaf_joint, linewidth=.5, alpha=1, color='k')
	f2ax4.axvline(x=0)
 	f2ax4.set_xlabel('$\Delta a_f/\\bar{a_f}$', fontsize='large', labelpad=10)
 	f2ax4.set_ylabel('$P(\Delta a_f/\\bar{a_f})$', fontsize='large', labelpad=10)
 	f2ax4.set_xlim(-0.5, 0.5)
 	fig2.tight_layout()
 	fig2.savefig('%s/P_epsilon_sigma_1D_cat%d.pdf' %(out_dir, cat_id))
 	fig2.savefig('%s/P_epsilon_sigma_1D_cat%d.png' %(out_dir, cat_id), dpi=200)

	# plot the 1-sigma and 2-sigma errors against the injections
 	plt.figure(figsize=(7,4.5))
 	ax = plt.subplot(211)
 	ax.set_xticklabels([])
 	plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left1_dMfbyMf_joint, right1_dMfbyMf_joint, color='#cc0000', alpha=0.8) #, label='68\%')
	plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left2_dMfbyMf_joint, right2_dMfbyMf_joint, color='#cc0000', alpha=0.2) #, label='95\%')
 	plt.plot(1+np.arange(len(mean_dMfbyMf)), (mean_dMfbyMf), '.', marker='o', color='#ff7c4c', mew=0, ms=8)
 	plt.errorbar(1+np.arange(len(mean_dMfbyMf)), mean_dMfbyMf, yerr=(upper_error_dMfbyMf, lower_error_dMfbyMf), linestyle='None', alpha=1, color='#ff7c4c', lw=0.1, capsize=0)
 	plt.axhline(y=0., color='k', ls='--', lw=0.5)
 	plt.ylim(-0.4,0.4)
 	plt.xlim(1, MAX_NUM_EVENTS)
 	plt.ylabel('$\Delta M_f/\\bar{M_f}$', fontsize='xx-large',  labelpad=10)
 	plt.yticks(fontsize = 'large')
 	plt.subplot(212)
 	plt.fill_between(1+np.arange(len(mean_dafbyaf_joint)), left1_dafbyaf_joint, right1_dafbyaf_joint, color='#cc0000', alpha=0.8) #, label='68\%')
 	plt.fill_between(1+np.arange(len(mean_dafbyaf_joint)), left2_dafbyaf_joint, right2_dafbyaf_joint, color='#cc0000', alpha=0.2) #, label='95\%')
 	plt.plot(1+np.arange(len(mean_dafbyaf)), (mean_dafbyaf), '.', marker='o', color='#ff7c4c', mew=0, ms=8)
 	plt.errorbar(1+np.arange(len(mean_dafbyaf)), mean_dafbyaf, yerr=(upper_error_dafbyaf, lower_error_dafbyaf), linestyle='None', alpha=1, color='#ff7c4c', lw=0.1, capsize=0)
 	plt.axhline(y=0., color='k', ls='--', lw=0.5)
 	plt.ylim(-0.4,0.4)
 	plt.xlim(1, MAX_NUM_EVENTS)
 	plt.xticks(fontsize = 'large')
	plt.yticks(fontsize = 'large')
 	plt.xlabel('Number of events', fontsize='xx-large',  labelpad=10)
 	plt.ylabel('$\Delta a_f/\\bar{a_f}$', fontsize='xx-large',  labelpad=10)
 	plt.tight_layout()
 	plt.savefig('%s/P_epsilon_sigma_marg1D_cat%d.pdf' %(out_dir,cat_id))
 	plt.savefig('%s/P_epsilon_sigma_marg1D_cat%d.png' %(out_dir, cat_id), dpi=200)

	plt.close()

  # make P-P plots
  # compute the histograms of the credible intervals
  plt.figure(figsize=(9.5,3.5))
  Nbins = 100
  N_cat = 25
  edges = np.linspace(0, 1, Nbins)

  # make p-p plot from different subset of samples
  for cat_id in range(N_cat):
  	  np.random.shuffle(gr_conf_arr)
  	  np.random.shuffle(gr_conf_eps_arr)
          np.random.shuffle(gr_conf_sigma_arr)

          if cat_id == N_cat-1:
                alph, col, N_sampl, linw = 1., 'r', len(gr_conf_arr), 2
          else:
                alph, col, N_sampl, linw = 0.2, 'k', 50, 1


	  ###FIXME####
	  gr_conf_arr[np.isnan(gr_conf_arr)] = 0
	  gr_conf_eps_arr[np.isnan(gr_conf_eps_arr)] = 0
	  gr_conf_sigma_arr[np.isnan(gr_conf_sigma_arr)] = 0
	  ###FIXME####

          Hepssig, edges_es = np.histogram(gr_conf_arr[:N_sampl], bins=edges, normed=True)
	  Heps, edges_eps = np.histogram(gr_conf_eps_arr[:N_sampl], bins=edges, normed=True)
          Hsig, edges_sig = np.histogram(gr_conf_sigma_arr[:N_sampl], bins=edges, normed=True)

          # compute cumulatives
          cdf_es = np.cumsum(Hepssig)*(edges_es[1]-edges_es[0])
          cdf_eps = np.cumsum(Heps)*(edges_eps[1]-edges_eps[0])
          cdf_sig = np.cumsum(Hsig)*(edges_sig[1]-edges_sig[0])

          # plot them
          plt.subplot(131)
          plt.plot(edges_es[1:], cdf_es, color=col, alpha=alph, lw=linw)
          plt.subplot(132)
          plt.plot(edges_eps[1:], cdf_eps, color=col, alpha=alph, lw=linw)
          plt.subplot(133)
          plt.plot(edges_sig[1:], cdf_sig, color=col, alpha=alph, lw=linw)
  plt.subplot(131)
  plt.plot(edges_es[1:],edges_es[1:], 'k', ls='--', lw=1)
  plt.xlabel('GR credible level', fontsize='large',  labelpad=10)
  plt.ylabel('Fraction of events (cumulative)', fontsize='large',  labelpad=10)
  plt.xlim(0,1)
  plt.ylim(0,1)
  plt.title('$(\Delta M_f/\\bar{M}_f, \Delta a_f/\\bar{a}_f)$')
  plt.subplot(132)
  plt.plot(edges_es[1:],edges_es[1:], 'k', ls='--', lw=1)
  plt.xlabel('GR credible level', fontsize='large',  labelpad=10)
  plt.xlim(0,1)
  plt.ylim(0,1)
  plt.title('$\Delta M_f/\\bar{M}_f$')
  plt.subplot(133)
  plt.plot(edges_es[1:],edges_es[1:], 'k', ls='--', lw=1)
  plt.xlabel('GR credible level', fontsize='large',  labelpad=10)
  plt.title('$\Delta a_f/\\bar{a}_f$')
  plt.xlim(0,1)
  plt.ylim(0,1)
  plt.tight_layout()
  plt.savefig('%s/PP_plots.png' %(out_dir), dpi=200)
  plt.savefig('%s/PP_plots.pdf' %(out_dir))
  plt.close()


if __name__ == '__main__':
	post_loc_list = np.genfromtxt('post_loc_list.txt', unpack=True, dtype='string')
	out_dir = '/home/rahul.kashyap/testGR_IR/runs/sim_injections_NOKDE/comb_Posterior2'
	prior_file = 'Prior_normbymean_nonprecspin_Healy2014_m12min3.0_m12max300.00_a12zmin0.00_a12zmax0.99_epsilon1.0_sigma1.0_401bins_m1m2a1a2'
	snr_thresh, GR_conf_thresh, q_thresh, mtot_thresh, Nbins, MAX_NUM_EVENTS = 8., 100., 20., 150., 401, 50

	calc_joint_posterior_module(post_loc_list, out_dir, prior_file, snr_thresh, GR_conf_thresh, q_thresh, mtot_thresh, Nbins, MAX_NUM_EVENTS)
