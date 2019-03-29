"""

Compute the joint posterior on Delta Mf/Mf and Delta chif/chif from multiple simulated BBH events

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

inj_no_list = np.arange(1, 5000, 1)

SNR_THRESH = 8.
GR_CONF_THRESH = 100.
Q_THRESH = 20.
MTOT_THRESH = 150.
Nbins = 401
MAX_NUM_EVENTS = 4.
SPIN_THRESH = 1. 

norm_by = 'norm_by_mean'
tag = 'snrthresh%d_GRthresh_%d_Mthresh_%d_qthresh_%d_chithresh%3.2f_4_events_snr_around_40' %(SNR_THRESH, GR_CONF_THRESH, MTOT_THRESH, Q_THRESH, SPIN_THRESH)

inj_type = 'GR'
#inj_type = 'modGR_a2_20'
#inj_type = 'IHES_GR'


if inj_type == 'GR':

	# using the same defn as above; now usign the same prior as lalinference; ie, uniform in component masses and spins 
	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-03-31/altdef/Nbins401'	
	out_dir = './CQG_response/%s'%tag

elif inj_type == 'modGR_a2_20':
	post_loc = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/modGR_a2_20/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-30_3det/uniform_compmass_spins_comoving_volume/imrtgr/2017-03-25/%s/Nbins%d/imrtgr' %(norm_by, Nbins)
	post_loc = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/modGR_a2_20/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-30_3det/uniform_compmass_spins_comoving_volume/imrtgr/2017-03-27/%s/Nbins%d/imrtgr' %(norm_by, Nbins)
	post_loc = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/modGR_a2_20/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-30_3det/uniform_compmass_spins_comoving_volume/imrtgr/2017-03-27/norm_by_mean/gfsmoothedprior_maxspin0.98_Nbins401/imrtgr'

	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/modGR_a2_20/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-03-31/altdef/Nbins401' # using the same defn as above; now usign the same prior as lalinference; ie, uniform in component masses and spins 
	out_dir  = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/modGR_a2_20/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-30_3det/uniform_compmass_spins_comoving_volume/imrtgr/2017-03-31/normbymean_gaussiansmoothed401bins_range2_noMfafprior/Nbins%d/%s/' %(Nbins, tag)

elif inj_type == 'IHES_GR':
	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/GR_20170407/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-03-31/altdef/Nbins401'
	out_dir = '/home/ajith/working/cbc/testGR_IR/runs/simulations/gaussian_noise/IHES_GR/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/20170407/imrtestgr_nonprecspin_Healy2014_normbymean_gaussiansmoothed401bins_range2_final_noMfafprior/imrtgr/2017-04-15/Nbins%d/%s/' %(Nbins, tag)

#prior_file = '../data/SemiAnalyticPrior_dMfbyMf_dafbyaf_Mfmin_1.0_Mfmax_300.0_afmin0_afmax1_dMfbyMflim2.0_dafbyaflim2.0'
#prior_file = '/home/abhirup/Documents/Work/testGR_IR/data/Prior_normbymean_updateddef_Mfmin_1.0_Mfmax_500.0_chifmin0.0_chifmax1.0_dMfbyMflim2.0_dchifbychiflim2.0_401bins_eq54'
prior_file = '/home/abhirup/Documents/Work/testGR_IR/data/Prior_normbymean_eq54_nonprecspin_Healy2014_m12min1.0_m12max300.00_a12zmin-0.98_a12zmax0.98_epsilon2.0_sigma2.0_aligned_401bins_m1m2a1za2z_gaussiansmoothed'

inj_file = '/home/abhirup/Documents/Work/testGR_IR/injections/popsynth_injections_SEOBNRv2_ROM_DoubleSpinthreePointFivePN.txt'
m1_inj_list, m2_inj_list, chi1_inj_list, chi2_inj_list = np.loadtxt(inj_file, usecols=(0,1,6,7), unpack=True)

mark_size = 8

# create epsilon and sigma vectors. The posteriors will be interpolated to these points 
epsilon_vec = np.linspace(-2, 2, Nbins)
sigma_vec = np.linspace(-2, 2, Nbins)

# create outiput directory 
os.system('mkdir -p %s'%out_dir) 
os.system('cp %s %s' %(__file__, out_dir))

# read the prior file and evaluate the prior 
f = gzip.open(prior_file+".pklz",'rb')
P_dMfbyMf_dchifbychif_pr_interp_obj = pickle.load(f)

for cat_id in range(3):
 print '# cat_id = %d' %cat_id 
 num_successful_events = 1

 # initialize variables 
 P_dMfbyMf_dchifbychif_joint = np.ones((Nbins, Nbins))
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
 gr_conf_eps_arr = np.array([])
 gr_conf_sigma_arr = np.array([])
 inj_no_arr = np.array([])
 Mf_inj_arr = np.array([])
 chif_inj_arr = np.array([])
 m1_inj_arr = np.array([])
 m2_inj_arr = np.array([])
 chi1_inj_arr = np.array([])
 chi2_inj_arr = np.array([])

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
 np.random.shuffle(inj_no_list) 

 for inj_no in inj_no_list:

	print '... injection: %d' %inj_no 

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

	# posterior file names 
	if norm_by == 'norm_by_insp': 
		post_file = '%s/data/P_MfRbyMfI_chifRbychifI.dat.gz' %in_dir
		epsilon_bins_file = '%s/data/MfRbyMfI_vec.dat.gz' %in_dir
		sigma_bins_file = '%s/data/chifRbychifI_vec.dat.gz' %in_dir
		gr_conf_file = '%s/data/GR_confidence_inspnorm.txt' %in_dir
	elif norm_by == 'norm_by_imr':
		post_file = '%s/data/P_dMfbyMf_dchifbychif.dat.gz' %in_dir
		epsilon_bins_file = '%s/data/dMfbyMf_vec.dat.gz' %in_dir
		sigma_bins_file = '%s/data/dchifbychif_vec.dat.gz' %in_dir
		gr_conf_file = '%s/data/GR_confidence.txt' %in_dir
	elif norm_by == 'norm_by_mean':
		#post_file = '%s/data/P_dMfbyMfmean_dchifbychifmean.dat.gz' %in_dir
		post_file = '%s/data/P_dMfbyMf_dchifbychif.dat.gz' %in_dir
		epsilon_bins_file = '%s/data/dMfbyMf_vec.dat.gz' %in_dir
		sigma_bins_file = '%s/data/dchifbychif_vec.dat.gz' %in_dir
		#gr_conf_file = '%s/data/GR_confidence_normbymean.txt' %in_dir
		gr_conf_file = '%s/data/GR_confidence.txt' %in_dir

	if os.path.isfile(post_file) is not True or os.path.isfile(epsilon_bins_file) is not True:
		print '...... data not found'
	else: 
		print '...... found data'

		# read the text files containing the optimal snr 
		insp_summary_stats_loc = '%s/lalinf_insp/summary_statistics.dat'%in_dir
		ring_summary_stats_loc = '%s/lalinf_ring/summary_statistics.dat'%in_dir
		imr_summary_stats_loc = '%s/lalinf_imr/summary_statistics.dat'%in_dir

		# read the optimal from each of the runs 
		insp_snr, insp_h1_snr, insp_l1_snr, insp_v1_snr, insp_m1_rec, insp_m2_rec, insp_chi1_rec, insp_chi2_rec, insp_Mf_rec, insp_chif_rec = optimal_snr_module(insp_summary_stats_loc)
		ring_snr, ring_h1_snr, ring_l1_snr, ring_v1_snr, ring_m1_rec, ring_m2_rec, ring_chi1_rec, ring_chi2_rec, ring_Mf_rec, ring_chif_rec = optimal_snr_module(ring_summary_stats_loc)
		imr_snr, imr_h1_snr, imr_l1_snr, imr_v1_snr, imr_m1_rec, imr_m2_rec, imr_chi1_rec, imr_chi2_rec, imr_Mf_rec, imr_chif_rec = optimal_snr_module(imr_summary_stats_loc)

		# read gr confidence value
		gr_conf_level = np.loadtxt(gr_conf_file)

		# apply some threshold 
		if insp_snr < SNR_THRESH or ring_snr < SNR_THRESH or imr_snr < 35. or imr_snr > 45. or q_inj > Q_THRESH or M_inj >  MTOT_THRESH or gr_conf_level > GR_CONF_THRESH/100. or abs(chi1_inj) > SPIN_THRESH or abs(chi2_inj) > SPIN_THRESH:
			print '...... did not pass thresholds'
		else: 
			print '...... passed thresholds'

			inj_no_arr = np.append(inj_no_arr, inj_no)
			Mf_inj_arr = np.append(Mf_inj_arr, Mf_inj)
			chif_inj_arr = np.append(chif_inj_arr, chif_inj)
			m1_inj_arr = np.append(m1_inj_arr, m1_inj)
			m2_inj_arr = np.append(m2_inj_arr, m2_inj)
			chi1_inj_arr = np.append(chi1_inj_arr, chi1_inj)
			chi2_inj_arr = np.append(chi2_inj_arr, chi2_inj)

			# read the posterior data for this event 
			P_dMfbyMf_dchifbychif_tmp = np.loadtxt(post_file)
			epsilon_vec_tmp = np.loadtxt(epsilon_bins_file)
			sigma_vec_tmp = np.loadtxt(sigma_bins_file)

			# interpolate the posterior to a fixed grid 
			P_dMfbyMf_dchifbychif_obj = scipy.interpolate.interp2d(epsilon_vec_tmp, sigma_vec_tmp, P_dMfbyMf_dchifbychif_tmp, fill_value=0., bounds_error=False)
			P_dMfbyMf_dchifbychif = P_dMfbyMf_dchifbychif_obj(epsilon_vec, sigma_vec)

			print '...... read posterior files'

			# if the posteriors are for (Mf_R/Mf_I, af_R/af_I) convert them to (1-Mf_R/Mf_I, 1-af_R/af_I)
			if norm_by == 'norm_by_insp': 
				epsilon_vec = 1.-epsilon_vec
				sigma_vec = 1.-sigma_vec

			# evaluate the prior at these bins 
			P_dMfbyMf_dchifbychif_pr = P_dMfbyMf_dchifbychif_pr_interp_obj(epsilon_vec, sigma_vec)	

			# Calculation of confidence levels in the 2D posterior 
			s1_v1v2 = ba.nsigma_value(P_dMfbyMf_dchifbychif, 0.68)
			s2_v1v2 = ba.nsigma_value(P_dMfbyMf_dchifbychif, 0.95)

			# plot the 2D 1-sigma contour from the current posteriors 
			if num_successful_events < MAX_NUM_EVENTS:
				ax1.contour(epsilon_vec,sigma_vec, gf(P_dMfbyMf_dchifbychif), levels=(s1_v1v2,), linewidths=(.5,), alpha=0.5, colors='#ff7c4c')

			# compute the confidence region corresponding to the GR value (delta_Mf/Mf = 0, delta_chif/chif = 0). 
			# the 'confidence' class is defined on top of this script. This is computed using a flat prior in delta_Mf/Mf, delta_chif/chif
			P_dMfbyMf_dchifbychif_flatprior = P_dMfbyMf_dchifbychif/P_dMfbyMf_dchifbychif_pr
			#conf_v1v2 = confidence(P_dMfbyMf_dchifbychif_flatprior)
			conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
			gr_height = P_dMfbyMf_dchifbychif[np.argmin(abs(epsilon_vec)), np.argmin(abs(sigma_vec))] # taking value closest to (0,0)
			gr_conf_level = conf_v1v2.level_from_height(gr_height)
			gr_conf_arr = np.append(gr_conf_arr, gr_conf_level)

			# removing nans and infinities 
			P_dMfbyMf_dchifbychif[np.isnan(P_dMfbyMf_dchifbychif)] = 0.
			P_dMfbyMf_dchifbychif[np.isinf(P_dMfbyMf_dchifbychif)] = 0.

			# replace all zeros in the posterior vector by a small number. 
			# This is to avoid multiplication by zeros while combining the posteiror 
			zidx = np.where(P_dMfbyMf_dchifbychif == 0.) 
			P_dMfbyMf_dchifbychif[zidx] = 1e-5 

			# Joint Probability distribution computation -- from second event onwards, divide the new posterior by the prior so that we are multiplying the 
			# the previous posterior by the current likelihood 
			if num_successful_events == 1: 
				P_dMfbyMf_dchifbychif_joint = P_dMfbyMf_dchifbychif
			elif num_successful_events > 1: 
				P_dMfbyMf_dchifbychif_joint = P_dMfbyMf_dchifbychif_joint*P_dMfbyMf_dchifbychif/P_dMfbyMf_dchifbychif_pr
			else: 
				raise ValueError('num_successful_events seems to be less than one!')

			dx = np.mean(np.diff(epsilon_vec))
			dy = np.mean(np.diff(sigma_vec))
			P_dMfbyMf_dchifbychif_joint /= np.sum(P_dMfbyMf_dchifbychif_joint) * dx * dy	# normalization 
				
			print '...... computed combined posteriors' 

			# Marginalization to one-dimensional joint_posteriors
			P_dMfbyMf_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=0) * dy
			P_dchifbychif_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=1) * dx
			P_dMfbyMf = np.sum(P_dMfbyMf_dchifbychif, axis=0) * dy
			P_dchifbychif = np.sum(P_dMfbyMf_dchifbychif, axis=1) * dx

			# calculate hte credible levels in the 1D posterior 
			conf_eps = confidence(P_dMfbyMf)
			gr_height_eps = P_dMfbyMf[np.argmin(abs(epsilon_vec))]
			gr_conf_level_eps = conf_eps.level_from_height(gr_height_eps)
			gr_conf_eps_arr = np.append(gr_conf_eps_arr, gr_conf_level_eps)

			conf_sigma = confidence(P_dchifbychif)
			gr_height_sigma = P_dchifbychif[np.argmin(abs(sigma_vec))]
			gr_conf_level_sigma = conf_sigma.level_from_height(gr_height_sigma)
			gr_conf_sigma_arr = np.append(gr_conf_sigma_arr, gr_conf_level_sigma)

			# plot the 1D marginalized posterior on dMfbyMf and dafbyaf
 			f2ax1.plot(epsilon_vec, P_dMfbyMf, linewidth=.5, alpha=0.5, color='#ff7c4c')
 			f2ax2.plot(sigma_vec, P_dchifbychif, linewidth=.5, alpha=0.5, color='#ff7c4c')

			# Calculation of credible intervals in the 2D posterior 
			s1_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint, 0.68)
			s2_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint, 0.95)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (joint) 
			s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_conf_intervals_in_1d(P_dMfbyMf_joint, epsilon_vec)
			s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_conf_intervals_in_1d(P_dchifbychif_joint, sigma_vec)

			# calcualte the confidence levels and intervals in the marginalized 1d posteriors (current injection) 
			s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_conf_intervals_in_1d(P_dMfbyMf, epsilon_vec)
			s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_conf_intervals_in_1d(P_dchifbychif, sigma_vec)

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
			mean_dMfbyMf_joint = np.append(mean_dMfbyMf_joint, np.average(epsilon_vec, weights=P_dMfbyMf_joint))
			mean_dchifbychif_joint = np.append(mean_dchifbychif_joint, np.average(sigma_vec, weights=P_dchifbychif_joint))

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
			mean_dMfbyMf = np.append(mean_dMfbyMf, np.average(epsilon_vec, weights=P_dMfbyMf))
			mean_dchifbychif = np.append(mean_dchifbychif, np.average(sigma_vec, weights=P_dchifbychif))

			# plot the combined posteriors from 10, 25, 50 events
			if num_successful_events == int(MAX_NUM_EVENTS/2) or num_successful_events == int(MAX_NUM_EVENTS/5) or num_successful_events == int(MAX_NUM_EVENTS):

				ax1.contour(epsilon_vec,sigma_vec, P_dMfbyMf_dchifbychif_joint, levels=(s1_v1v2_joint,), linewidths=(1.5,), alpha=num_successful_events/MAX_NUM_EVENTS, colors='#cc0000')
				ax2.contour(epsilon_vec,sigma_vec, P_dMfbyMf_dchifbychif_joint, levels=(s1_v1v2_joint,), linewidths=(1.5,), alpha=num_successful_events/MAX_NUM_EVENTS, colors='#cc0000')

				# plot the 1D marginalized posterior on dMfbyMf and dafbyaf
 				f2ax3.plot(epsilon_vec, P_dMfbyMf_joint, linewidth=1.5, alpha=num_successful_events/MAX_NUM_EVENTS, color='#cc0000')
 				f2ax4.plot(sigma_vec, P_dchifbychif_joint, linewidth=1.5, alpha=num_successful_events/MAX_NUM_EVENTS, color='#cc0000')

			num_successful_events += 1 

 #########################################################################
 ######## plot the confidence intervals against the event number  ########
 #########################################################################
 
 upper_error_dMfbyMf = mean_dMfbyMf - left1_dMfbyMf
 lower_error_dMfbyMf = right1_dMfbyMf - mean_dMfbyMf
 upper_error_dchifbychif = mean_dchifbychif - left1_dchifbychif
 lower_error_dchifbychif = right1_dchifbychif - mean_dchifbychif  
 
 # plot 2D posteriors 
 ax1.plot(0, 0, 'k+', mew=1.5, ms=12, alpha=1)
 ax1.set_xlabel('$\Delta M_f/\\bar{M_f}$', fontsize='xx-large', labelpad=10)
 ax1.set_ylabel('$\Delta a_f/\\bar{a_f}$', fontsize='xx-large', labelpad=10)
 ax1.set_xticks([-2, -1, 0, 1, 2])
 ax1.set_yticks([-2, -1, 0, 1, 2])
 
 ax2.plot(0, 0, 'k+', mew=1.5, ms=12, alpha=1)
 if inj_type == 'GR':
 	ax2.set_xlim([-0.1, 0.1]); ax2.set_ylim([-0.1, 0.1])
 else: 
 	ax2.set_xlim([-0.2, 0.2]); ax2.set_ylim([-0.2, 0.2])
 	ax2.set_xticks([-.2, 0, .2])
 	ax2.set_yticks([-.2, 0, .2])
 fig1.tight_layout()
 fig1.savefig('%s/P_epsilon_sigma_2D_%s_cat%d.pdf' %(out_dir, inj_type, cat_id))
 fig1.savefig('%s/P_epsilon_sigma_2D_%s_cat%d.png' %(out_dir, inj_type, cat_id), dpi=200)
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
 f2ax4.plot(sigma_vec, P_dchifbychif_joint, linewidth=.5, alpha=1, color='k')
 f2ax4.axvline(x=0)
 f2ax4.set_xlabel('$\Delta a_f/\\bar{a_f}$', fontsize='large', labelpad=10)
 f2ax4.set_ylabel('$P(\Delta a_f/\\bar{a_f})$', fontsize='large', labelpad=10)
 f2ax4.set_xlim(-0.5, 0.5)
 fig2.tight_layout()
 fig2.savefig('%s/P_epsilon_sigma_1D_%s_cat%d.pdf' %(out_dir, inj_type, cat_id))
 fig2.savefig('%s/P_epsilon_sigma_1D_%s_cat%d.png' %(out_dir, inj_type, cat_id), dpi=200)

 # plot the GR credible levels 
 plt.figure(figsize=(8,5))
 plt.plot(1+np.arange(len(gr_conf_arr)), gr_conf_arr*100, 'go-', label='GR cred level ($\epsilon, \sigma$)', ms=mark_size, lw=1, mew=0)
 plt.plot(1+np.arange(len(gr_conf_arr)), gr_conf_eps_arr*100, 'ro-', label='GR cred level ($\epsilon$)', ms=mark_size, lw=1, mew=0)
 plt.plot(1+np.arange(len(gr_conf_arr)), gr_conf_sigma_arr*100, 'bo-', label='GR cred level ($\sigma$)', ms=mark_size, lw=1, mew=0)
 plt.legend(loc='lower right', fontsize=5, frameon=False)
 plt.xlim(1, len(mean_dMfbyMf))
 plt.ylim(0,100)
 plt.ylabel('mean $\\rho_\mathrm{opt}$')
 plt.grid()
 plt.tight_layout()
 plt.savefig('%s/GR_confidence_levels_catid%d.png'%(out_dir, cat_id))
 plt.close()
 
 # plot the 1-sigma and 2-sigma errors against the injections 
 plt.figure(figsize=(7,4.5))
 ax = plt.subplot(211)
 ax.set_xticklabels([])
 plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left1_dMfbyMf_joint, right1_dMfbyMf_joint, color='#cc0000', alpha=0.8) #, label='68\%')
 plt.fill_between(1+np.arange(len(mean_dMfbyMf_joint)), left2_dMfbyMf_joint, right2_dMfbyMf_joint, color='#cc0000', alpha=0.2) #, label='95\%')
 plt.plot(1+np.arange(len(mean_dMfbyMf)), (mean_dMfbyMf), '.', marker='o', color='#ff7c4c', mew=0, ms=mark_size)
 plt.errorbar(1+np.arange(len(mean_dMfbyMf)), mean_dMfbyMf, yerr=(upper_error_dMfbyMf, lower_error_dMfbyMf), linestyle='None', alpha=1, color='#ff7c4c', lw=0.1, capsize=0)
 plt.axhline(y=0., color='k', ls='--', lw=0.5)
 plt.ylim(-0.4,0.4)
 plt.xlim(1, MAX_NUM_EVENTS)
 plt.ylabel('$\Delta M_f/\\bar{M_f}$', fontsize='xx-large',  labelpad=10)
 plt.yticks(fontsize = 'large')
 plt.subplot(212)
 plt.fill_between(1+np.arange(len(mean_dchifbychif_joint)), left1_dchifbychif_joint, right1_dchifbychif_joint, color='#cc0000', alpha=0.8) #, label='68\%')
 plt.fill_between(1+np.arange(len(mean_dchifbychif_joint)), left2_dchifbychif_joint, right2_dchifbychif_joint, color='#cc0000', alpha=0.2) #, label='95\%')
 plt.plot(1+np.arange(len(mean_dchifbychif)), (mean_dchifbychif), '.', marker='o', color='#ff7c4c', mew=0, ms=mark_size)
 plt.errorbar(1+np.arange(len(mean_dchifbychif)), mean_dchifbychif, yerr=(upper_error_dchifbychif, lower_error_dchifbychif), linestyle='None', alpha=1, color='#ff7c4c', lw=0.1, capsize=0)
 plt.axhline(y=0., color='k', ls='--', lw=0.5)
 plt.ylim(-0.4,0.4)
 plt.xlim(1, MAX_NUM_EVENTS)
 plt.xticks(fontsize = 'large')
 plt.yticks(fontsize = 'large')
 plt.xlabel('Number of events', fontsize='xx-large',  labelpad=10)
 plt.ylabel('$\Delta a_f/\\bar{a_f}$', fontsize='xx-large',  labelpad=10)
 plt.tight_layout()
 plt.savefig('%s/P_epsilon_sigma_marg1D_%s_cat%d.pdf' %(out_dir, inj_type, cat_id))
 plt.savefig('%s/P_epsilon_sigma_marg1D_%s_cat%d.png' %(out_dir, inj_type, cat_id), dpi=200)
 plt.close()

# statistical errors 
plt.figure(figsize=(5,5))
inj_vec = 1+np.arange(len(mean_dMfbyMf_joint))
plt.loglog(inj_vec, abs(left1_dMfbyMf_joint-right1_dMfbyMf_joint), 'k', label='$\Delta M_f/M_f$')
plt.loglog(inj_vec, abs(left1_dchifbychif_joint-right1_dchifbychif_joint), 'orange', label='$\Delta \chi_f/\chi_f$')
plt.ylabel('Combined 68\% confidence')
plt.xlabel('Number of events')
plt.ylim(1e-2,1)
plt.tight_layout()
plt.savefig('%s/statistical_errors.png'%out_dir, dpi=300)
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
plt.savefig('%s/PP_plots_%s.png' %(out_dir, inj_type), dpi=200)
plt.savefig('%s/PP_plots_%s.pdf' %(out_dir, inj_type))
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
plt.savefig('%s/mean_and_error_vs_numinj.png' %out_dir, dpi=300)
plt.close()

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
plt.savefig('%s/error_vs_numinj.pdf' %out_dir)
plt.close()


#########################################################################
######## save relevant data  ########
#########################################################################

np.savetxt('%s/gr_credible_levels.dat' %out_dir, np.column_stack([inj_no_arr, gr_conf_arr, gr_conf_eps_arr, gr_conf_sigma_arr]))
