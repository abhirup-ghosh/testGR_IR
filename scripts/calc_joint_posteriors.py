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

def optimal_snr_module(filename):
	data = np.genfromtxt(filename, dtype=None, names=True, usecols=(0,1,2,3,4,5,6))
	var_names = [d[0] for d in data]
	stat_names = data.dtype.names
	optimal_snr = data[var_names.index('optimal_snr')][stat_names.index('mean')+1]	
	h1_optimal_snr = data[var_names.index('h1_optimal_snr')][stat_names.index('mean')+1]
	l1_optimal_snr = data[var_names.index('l1_optimal_snr')][stat_names.index('mean')+1]
	v1_optimal_snr = data[var_names.index('v1_optimal_snr')][stat_names.index('mean')+1]
	m1_rec = data[var_names.index('m1')][stat_names.index('mean')+1]
	m2_rec = data[var_names.index('m2')][stat_names.index('mean')+1]
	a1z_rec = data[var_names.index('a1z')][stat_names.index('mean')+1]
	a2z_rec = data[var_names.index('a2z')][stat_names.index('mean')+1]
	Mf_rec = data[var_names.index('mf')][stat_names.index('mean')+1]
	af_rec = data[var_names.index('af')][stat_names.index('mean')+1]
	return optimal_snr, h1_optimal_snr, l1_optimal_snr, v1_optimal_snr, m1_rec, m2_rec, a1z_rec, a2z_rec, Mf_rec, af_rec

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

inj_no_list = np.linspace(1,5000,5000)

a_vec = np.linspace(0.5, 1, len(inj_no_list))
cat_id = 0 
#random.shuffle(inj_no_list)
SNR_THRESH = 8.
GR_CONF_THRESH = 100.
Q_THRESH = 20.
MTOT_THRESH = 200.
prior_correction = True
	
#inj_dir_root = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/modGR_a2_20'
inj_dir_root = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume'
#in_dir_root = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_-2.0_2.0_w_priorcorr'%inj_dir_root
in_dir_root = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_bins_-2.0_2.0_w_priorcorr'%inj_dir_root
inj_file = '../injections/popsynth_injections_SEOBNRv2_ROM_DoubleSpinthreePointFivePN.txt'
m1_inj_list, m2_inj_list, chi1_inj_list, chi2_inj_list = np.loadtxt(inj_file, usecols=(0,1,6,7), unpack=True)

out_dir = '%s/snrthresh%d_GRthresh_%d_Mthresh_%d_qthresh_%d_catalog%d_uniformdMfbyMfdchifbychifprior' %(in_dir_root, SNR_THRESH, GR_CONF_THRESH, MTOT_THRESH, Q_THRESH, cat_id)
mark_size = 2

# create output directory 
os.system('mkdir -p %s'%out_dir) 

# initialize variables 
N_bins = 401
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
insp_h1_snr_arr = np.array([])
insp_l1_snr_arr = np.array([])
insp_v1_snr_arr = np.array([])
ring_snr_arr = np.array([])
ring_h1_snr_arr = np.array([])
ring_l1_snr_arr = np.array([])
ring_v1_snr_arr = np.array([])
imr_snr_arr = np.array([])
imr_h1_snr_arr = np.array([])
imr_l1_snr_arr = np.array([])
imr_v1_snr_arr = np.array([])
insp_Mf_rec_arr = np.array([])
insp_chif_rec_arr = np.array([])
insp_m1_rec_arr = np.array([])
insp_m2_rec_arr = np.array([])
insp_chi1_rec_arr = np.array([])
insp_chi2_rec_arr = np.array([])
ring_Mf_rec_arr = np.array([])
ring_chif_rec_arr = np.array([])
ring_m1_rec_arr = np.array([])
ring_m2_rec_arr = np.array([])
ring_chi1_rec_arr = np.array([])
ring_chi2_rec_arr = np.array([])
imr_Mf_rec_arr = np.array([])
imr_chif_rec_arr = np.array([])
imr_m1_rec_arr = np.array([])
imr_m2_rec_arr = np.array([])
imr_chi1_rec_arr = np.array([])
imr_chi2_rec_arr = np.array([])
gr_conf_arr = np.array([])
inj_no_arr = np.array([])
Mf_inj_arr = np.array([])
chif_inj_arr = np.array([])
m1_inj_arr = np.array([])
m2_inj_arr = np.array([])
chi1_inj_arr = np.array([])
chi2_inj_arr = np.array([])

shuffle(inj_no_list)

fig1 = plt.figure(figsize=(5, 5))
ax1 = fig1.add_subplot(111)
for inj_no in inj_no_list:
        m1_inj, m2_inj, chi1_inj, chi2_inj = m1_inj_list[inj_no-1], m2_inj_list[inj_no-1], chi1_inj_list[inj_no-1], chi2_inj_list[inj_no-1]
        #m1_inj, m2_inj, chi1_inj, chi2_inj = m1_inj_list[inj_no-1], m2_inj_list[inj_no-1], 0., 0.

        if m1_inj > m2_inj:
        	m1_inj, m2_inj = m2_inj, m1_inj
		chi1_inj, chi2_inj = chi2_inj, chi1_inj

	#in_dir = '%s/IHES_%04d' %(in_dir_root, inj_no)
	in_dir = '%s/injection_%d' %(in_dir_root, inj_no)
	M_inj = m1_inj+m2_inj
	q_inj = m1_inj/m2_inj
	Mf_inj, chif_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, fit_formula)
		

	if os.path.isfile('%s/data/P_dMfbyMf_dchifbychif.dat.gz' %in_dir) == True:

	  try:
		print '... found data'

		# read the text files containing the optimal snr 
		insp_summary_stats_loc = '%s/lalinf_insp/summary_statistics.dat'%in_dir
		ring_summary_stats_loc = '%s/lalinf_ring/summary_statistics.dat'%in_dir
		imr_summary_stats_loc = '%s/lalinf_imr/summary_statistics.dat'%in_dir

		# read the optimal from each of the runs 
		insp_snr, insp_h1_snr, insp_l1_snr, insp_v1_snr, insp_m1_rec, insp_m2_rec, insp_chi1_rec, insp_chi2_rec, insp_Mf_rec, insp_chif_rec = optimal_snr_module(insp_summary_stats_loc)
		ring_snr, ring_h1_snr, ring_l1_snr, ring_v1_snr, ring_m1_rec, ring_m2_rec, ring_chi1_rec, ring_chi2_rec, ring_Mf_rec, ring_chif_rec = optimal_snr_module(ring_summary_stats_loc)
		imr_snr, imr_h1_snr, imr_l1_snr, imr_v1_snr, imr_m1_rec, imr_m2_rec, imr_chi1_rec, imr_chi2_rec, imr_Mf_rec, imr_chif_rec = optimal_snr_module(imr_summary_stats_loc)

		# read gr confidence value
		gr_conf_level = np.loadtxt('%s/data/GR_confidence.txt'%in_dir)

		# apply some threshold 
		if insp_snr > SNR_THRESH and ring_snr > SNR_THRESH and q_inj < Q_THRESH and M_inj < MTOT_THRESH and gr_conf_level < GR_CONF_THRESH/100.:
			
			inj_no_arr = np.append(inj_no_arr, inj_no)
			Mf_inj_arr = np.append(Mf_inj_arr, Mf_inj)
			chif_inj_arr = np.append(chif_inj_arr, chif_inj)
			m1_inj_arr = np.append(m1_inj_arr, m1_inj)
			m2_inj_arr = np.append(m2_inj_arr, m2_inj)
			chi1_inj_arr = np.append(chi1_inj_arr, chi1_inj)
			chi2_inj_arr = np.append(chi2_inj_arr, chi2_inj)

			# append to the SNR vector (SNRs for each event) 
			insp_snr_arr = np.append(insp_snr_arr, insp_snr)
			insp_h1_snr_arr = np.append(insp_h1_snr_arr, insp_h1_snr)
			insp_l1_snr_arr = np.append(insp_l1_snr_arr, insp_l1_snr)
			insp_v1_snr_arr = np.append(insp_v1_snr_arr, insp_v1_snr)
			insp_Mf_rec_arr = np.append(insp_Mf_rec_arr, insp_Mf_rec)
			insp_chif_rec_arr = np.append(insp_chif_rec_arr, insp_chif_rec)
			insp_m1_rec_arr = np.append(insp_m1_rec_arr, insp_m1_rec)
			insp_m2_rec_arr = np.append(insp_m2_rec_arr, insp_m2_rec)
			insp_chi1_rec_arr = np.append(insp_chi1_rec_arr, insp_chi1_rec)
			insp_chi2_rec_arr = np.append(insp_chi2_rec_arr, insp_chi2_rec)
			ring_snr_arr = np.append(ring_snr_arr, ring_snr)
			ring_h1_snr_arr = np.append(ring_h1_snr_arr, ring_h1_snr)
			ring_l1_snr_arr = np.append(ring_l1_snr_arr, ring_l1_snr)
			ring_v1_snr_arr = np.append(ring_v1_snr_arr, ring_v1_snr)
			ring_Mf_rec_arr = np.append(ring_Mf_rec_arr, ring_Mf_rec)
                        ring_chif_rec_arr = np.append(ring_chif_rec_arr, ring_chif_rec)
                        ring_m1_rec_arr = np.append(ring_m1_rec_arr, ring_m1_rec)
                        ring_m2_rec_arr = np.append(ring_m2_rec_arr, ring_m2_rec)
                        ring_chi1_rec_arr = np.append(ring_chi1_rec_arr, ring_chi1_rec)
                        ring_chi2_rec_arr = np.append(ring_chi2_rec_arr, ring_chi2_rec)
			imr_snr_arr = np.append(imr_snr_arr, imr_snr)
			imr_h1_snr_arr = np.append(imr_h1_snr_arr, imr_h1_snr)
			imr_l1_snr_arr = np.append(imr_l1_snr_arr, imr_l1_snr)
			imr_v1_snr_arr = np.append(imr_v1_snr_arr, imr_v1_snr)
			imr_Mf_rec_arr = np.append(imr_Mf_rec_arr, imr_Mf_rec)
                        imr_chif_rec_arr = np.append(imr_chif_rec_arr, imr_chif_rec)
                        imr_m1_rec_arr = np.append(imr_m1_rec_arr, imr_m1_rec)
                        imr_m2_rec_arr = np.append(imr_m2_rec_arr, imr_m2_rec)
                        imr_chi1_rec_arr = np.append(imr_chi1_rec_arr, imr_chi1_rec)
                        imr_chi2_rec_arr = np.append(imr_chi2_rec_arr, imr_chi2_rec)

			# read the posterior data for this event 
			P_dMfbyMf_dchifbychif = np.loadtxt('%s/data/P_dMfbyMf_dchifbychif.dat.gz'%in_dir)
			dMfbyMf_vec = np.loadtxt('%s/data/dMfbyMf_vec.dat.gz'%in_dir)
			dchifbychif_vec = np.loadtxt('%s/data/dchifbychif_vec.dat.gz'%in_dir)

			# Calculation of confidence levels in the 2D posterior 
			s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dchifbychif, 0.68)
			s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dchifbychif, 0.95)

			# plot the 2D 1-sigma contour from the current posteriors 
                        ax1.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif), levels=(s1_v1v2,), linewidths=(0.1,), alpha=0.5, colors='#ff7c4c')
                        ax1.hold(True)

			print 'prior_correction = ', prior_correction

			# for the second event onwards, take out the prior -- the prior is included only once -- for the first event 
			if prior_correction == True: 
				prior_file = '../data/Prior_dMfbyMf_dchifbychif_Mfmin_1.0_Mfmax_500.0_chifmin-1.0_chifmax1.0_dMfbyMflim2.0_dchifbychiflim2.0'
				f = gzip.open(prior_file+".pklz",'rb')
				P_dMfbyMf_dchifbychif_pr_interp_obj = pickle.load(f)
				P_dMfbyMf_dchifbychif_pr = P_dMfbyMf_dchifbychif_pr_interp_obj(dMfbyMf_vec, dchifbychif_vec)	
				P_dMfbyMf_dchifbychif /= P_dMfbyMf_dchifbychif_pr

			# for the second event onwards, take out the prior -- the prior is included only once -- for the first event 
			prior_correction = True

			#compute the confidence region corresponding to the GR value (delta_Mf/Mf = 0, delta_chif/chif = 0). 
                        # the 'confidence' class is defined on top of this script 
                        conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
                        gr_height = P_dMfbyMf_dchifbychif[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dchifbychif_vec))] # taking value closest to (0,0)
                        gr_conf_level = conf_v1v2.level_from_height(gr_height)
                        gr_conf_arr = np.append(gr_conf_arr, gr_conf_level)

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

			print len(inj_no_arr)
			if len(inj_no_arr) == 50:
				P_dMfbyMf_dchifbychif_joint_50 = P_dMfbyMf_dchifbychif_joint
				s1_v1v2_joint_50 = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint_50, 0.68)
                        	s2_v1v2_joint_50 = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint_50, 0.95)

			if len(inj_no_arr) == 100:
				break

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

	  except IOError:
			print 'data not found'

#########################################################################
######## plot the confidence intervals against the event number  ########
#########################################################################

upper_error_dMfbyMf = mean_dMfbyMf - left1_dMfbyMf
lower_error_dMfbyMf = right1_dMfbyMf - mean_dMfbyMf
upper_error_dchifbychif = mean_dchifbychif - left1_dchifbychif
lower_error_dchifbychif = right1_dchifbychif - mean_dchifbychif  

ax1.contour(dMfbyMf_vec,dchifbychif_vec,gf(P_dMfbyMf_dchifbychif_joint), levels=(s1_v1v2_joint,s2_v1v2_joint), linewidths=(0.5,0.5), colors='#cc0000')
ax1.plot(0, 0, 'k+', mew=0.5)
plt.xticks(fontsize = 'large')
plt.yticks(fontsize = 'large')
plt.xlabel('$\Delta M_f/M_f$', fontsize='xx-large')
plt.ylabel('$\Delta a_f/a_f$', fontsize='xx-large')
plt.tight_layout()
left, bottom, width, height = [0.32, 0.7, 0.2, 0.2]
ax2 = plt.axes([left, bottom, width, height])
ax2.contour(dMfbyMf_vec,dchifbychif_vec,gf(P_dMfbyMf_dchifbychif_joint_50), levels=(s2_v1v2_joint_50,s2_v1v2_joint_50), linewidths=(0.5, 0.5), colors='orange')
ax2.contour(dMfbyMf_vec,dchifbychif_vec,gf(P_dMfbyMf_dchifbychif_joint), levels=(s2_v1v2_joint, s2_v1v2_joint), linewidths=(1, 1), colors='#cc0000')
ax2.plot(0, 0, 'k+', mew=0.5)
plt.xlim([-0.05, 0.1]); plt.ylim([-0.05, 0.1])
#plt.xlim([-0.1, 0.05]); plt.ylim([-0.2, 0.05])
plt.setp(ax2, xlabel='$\Delta M_f/M_f$', ylabel='$\Delta a_f/a_f$', xticks=[-0.05, 0, 0.1], yticks=[-0.05, 0, 0.1])
#plt.setp(ax2, xlabel='$\Delta M_f/M_f$', ylabel='$\Delta a_f/a_f$', xticks=[-0.1, 0, 0.05], yticks=[-0.2, 0, 0.05])
#plt.savefig('%s/dMfbyMf_dchifbychif_jointposterior_vs_events_rightpanel.pdf'%out_dir, dpi=300)
#plt.savefig('/home/abhirup/Documents/Work/imrtgr_detailed_paper/fig/dMfbyMf_dchifbychif_jointposterior_vs_events_IHESmodGR_rightpanel_realisation2.pdf', dpi=300)
plt.savefig('/home/abhirup/Documents/Work/imrtgr_detailed_paper/fig/dMfbyMf_dchifbychif_jointposterior_vs_events_GR_injections_LALSimulation_rightpanel_realisation2.pdf', dpi=300)

plt.figure(figsize=(8,5))
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
plt.tight_layout()
#plt.savefig('%s/GR_confidence_levels.pdf'%out_dir, dpi=300)

plt.figure(figsize=(8,5))
ax = plt.subplot(211)
ax.set_xticklabels([])
plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left1_dMfbyMf_joint, right1_dMfbyMf_joint, color='#cc0000', alpha=0.8) #, label='68\%')
plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left2_dMfbyMf_joint, right2_dMfbyMf_joint, color='#cc0000', alpha=0.2) #, label='95\%')
plt.plot(1+np.arange(len(mean_dMfbyMf)), (mean_dMfbyMf), '.', marker='o', color='#ff7c4c', mew=0, ms=mark_size)
plt.errorbar(1+np.arange(len(mean_dMfbyMf)), mean_dMfbyMf, yerr=(upper_error_dMfbyMf, lower_error_dMfbyMf), linestyle='None', alpha=1, color='#ff7c4c', lw=0.1, capsize=0)
plt.axhline(y=0., color='k', ls='--', lw=0.5)
plt.ylim(-0.4,0.4)
plt.xlim(1, len(mean_dMfbyMf))
plt.ylabel('$\Delta M_f/M_f$', fontsize='xx-large')
plt.yticks(fontsize = 'large')
plt.subplot(212)
plt.fill_between(1+np.arange(len(mean_dchifbychif_joint)), left1_dchifbychif_joint, right1_dchifbychif_joint, color='#cc0000', alpha=0.8) #, label='68\%')
plt.fill_between(1+np.arange(len(mean_dchifbychif_joint)), left2_dchifbychif_joint, right2_dchifbychif_joint, color='#cc0000', alpha=0.2) #, label='95\%')
plt.plot(1+np.arange(len(mean_dchifbychif)), (mean_dchifbychif), '.', marker='o', color='#ff7c4c', mew=0, ms=mark_size)
plt.errorbar(1+np.arange(len(mean_dchifbychif)), mean_dchifbychif, yerr=(upper_error_dchifbychif, lower_error_dchifbychif), linestyle='None', alpha=1, color='#ff7c4c', lw=0.1, capsize=0)
plt.axhline(y=0., color='k', ls='--', lw=0.5)
plt.ylim(-0.4,0.4)
plt.xlim(1, len(mean_dMfbyMf))
plt.xticks(fontsize = 'large')
plt.yticks(fontsize = 'large')
plt.xlabel('Number of events', fontsize='xx-large')
plt.ylabel('$\Delta a_f/a_f$', fontsize='xx-large')
plt.tight_layout()
#plt.savefig('%s/dMfbyMf_dchifbychif_jointposterior_vs_events_leftpanel.pdf'%out_dir, dpi=300)
#plt.savefig('/home/abhirup/Documents/Work/imrtgr_detailed_paper/fig/dMfbyMf_dchifbychif_jointposterior_vs_events_IHESmodGR_leftpanel_realisation2.pdf', dpi=300)
plt.savefig('/home/abhirup/Documents/Work/imrtgr_detailed_paper/fig/dMfbyMf_dchifbychif_jointposterior_vs_events_GR_injections_LALSimulation_leftpanel_realisation2.pdf', dpi=300)

# plot the 1-sigma and 2-sigma errors against the injections 
plt.figure(figsize=(5,5))
inj_vec = 1+np.arange(len(mean_dMfbyMf_joint))
plt.loglog(inj_vec, abs(left1_dMfbyMf_joint-right1_dMfbyMf_joint), 'k', label='$\Delta M_f/M_f$')
plt.loglog(inj_vec, abs(left1_dchifbychif_joint-right1_dchifbychif_joint), 'orange', label='$\Delta \chi_f/\chi_f$')
plt.ylabel('Combined 68\% confidence')
plt.xlabel('Number of events')
plt.ylim(1e-2,1)
plt.tight_layout()
#plt.savefig('%s/statistical_errors.png'%out_dir, dpi=300)

# make the P-P plot
px_vec = np.linspace(0, 1, 100)
py_vec = np.zeros(100)
num_inj = float(len(gr_conf_arr))

for ix, px in enumerate(px_vec):
	py_vec[ix] = len(np.where(gr_conf_arr <= px)[0])/num_inj

realisations = 50
	
plt.figure(figsize=(5,5))
for idx in range(realisations):
	shuffle(gr_conf_arr)
	gr_conf_rand = gr_conf_arr[:100]
	px_vec_rand = np.linspace(0, 1, 100)
	py_vec_rand = np.zeros(100)
	num_inj_rand = float(len(gr_conf_rand))
	for ix, px in enumerate(px_vec_rand):
	        py_vec_rand[ix] = len(np.where(gr_conf_rand <= px)[0])/num_inj_rand	
	plt.plot(px_vec_rand, py_vec_rand, 'k', alpha=0.1, lw=1)
	plt.hold(True)
plt.plot(px_vec, px_vec, 'c', ls='--')
plt.plot(px_vec, py_vec, 'r')
plt.grid()
plt.xlabel('GR confidence threshold', fontsize='large')
plt.ylabel('fraction of events with GR confidence \nbelow the threshold', fontsize='large')
plt.xticks(fontsize = 'large')
plt.yticks(fontsize = 'large')
plt.tight_layout()
#plt.savefig('%s/pp_plots.pdf' %out_dir, dpi=300)
#plt.savefig('/home/abhirup/Documents/Work/testGR_IR/papers/imrtgr_followup_paper/fig/pp_plot_IHESmodGR.pdf', dpi=300)
#plt.savefig('/home/abhirup/Documents/Work/testGR_IR/papers/imrtgr_followup_paper/fig/pp_plot_GR_injections_LALSimulation.pdf', dpi=300)
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

#plt.savefig('%s/mean_and_error_vs_numinj.png' %out_dir, dpi=300)

# plot the 1-sigma and 2-sigma errors against the injections 
inj_vec = 1+np.arange(len(mean_dMfbyMf_joint))
plt.figure(figsize=(2.6,2.5))
plt.loglog(inj_vec, abs(left1_dMfbyMf_joint-right1_dMfbyMf_joint), 'k', label='$\Delta M_f/M_f$')
plt.loglog(inj_vec, abs(left1_dchifbychif_joint-right1_dchifbychif_joint), 'orange', label='$\Delta \chi_f/\chi_f$')
plt.ylabel('Combined 68\% confidence')
plt.xlabel('Number of events')
plt.ylim(1e-2,1)
plt.tight_layout()
plt.legend(frameon=False)
#plt.savefig('%s/error_vs_numinj.pdf' %out_dir)


#########################################################################
######## save relevant data  ########
#########################################################################

np.savetxt('/home/abhirup/Documents/Work/imrtgr_detailed_paper/data/inj_no_LALSimulation_realisation2.dat', np.c_[inj_no_arr])
np.savetxt('/home/abhirup/Documents/Work/imrtgr_detailed_paper/data/gr_conf_LALSimulation_realisation2.dat', np.c_[gr_conf_arr])
