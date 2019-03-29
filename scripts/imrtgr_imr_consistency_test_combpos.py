"""
Compute the joint posterior on Delta Mf/Mf and Delta chif/chif from multiple simulated BBH events

A. Ghosh, P. Ajith, 2015-11-27
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage.filters as filter
import time
import os
import os.path
import commands
import imrtestgr as tgr
import pickle, gzip
from scipy import interpolate
import tgrplotsettings
import optparse as op
from matplotlib import rc
import matplotlib
matplotlib.rc('text.latex', preamble = '\usepackage{txfonts}')

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='times')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=1)
rc("lines", linewidth=2)
rc('axes', labelsize=10)
rc("axes", linewidth=0.5)
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)
rc('legend', fontsize=10)
rc('xtick.major', pad=6)
rc('ytick.major', pad=6)
rc('xtick.minor', size=5)
rc('ytick.minor', size=5)

# plottting style 
def set_tick_sizes(ax, major, minor):
  for l in ax.get_xticklines() + ax.get_yticklines():
    l.set_markersize(major)
  for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
    tick.tick1line.set_markersize(minor)
    tick.tick2line.set_markersize(minor)
  ax.xaxis.LABELPAD=10.
  ax.xaxis.OFFSETTEXTPAD=10.

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

# gaussian filter with 1-sigma softening 
def gf(P):
	return filter.gaussian_filter(P, sigma=1.0)

# compute 1-sigma confidence intervals in 1D 
def calc_conf_intervals_in_1d(P, x):

	# find the value of P corresponding to 68% and 95% confidence heights 
	conf = confidence(P)
	P_s1 = conf.height_from_level(0.68)
	P_s2 = conf.height_from_level(0.68)

	# calculation of condifence edges (values of x corresponding to the height s1 on the two sides) 
	x_s1_l = min(x[np.where(P >= P_s1)[0]])
	x_s1_r = max(x[np.where(P >= P_s1)[0]])

	# calculation of condifence edges (values of x corresponding to the height s2 on the two sides) 
	x_s2_l = min(x[np.where(P >= P_s2)[0]])
	x_s2_r = max(x[np.where(P >= P_s2)[0]])

	return P_s1, P_s2, x_s1_l, x_s1_r, x_s2_l, x_s2_r

# main program
if __name__ == "__main__":

	parser = op.OptionParser()
	parser.add_option("-d", "--postloc-list", dest="post_loc_list", action="append", help="paths to the all the relevant data folders")
	parser.add_option("-o","--out-dir", dest="out_dir", help="output directory")
	(options, args) = parser.parse_args()
	post_loc_list = options.post_loc_list
	out_dir = options.out_dir

	os.system('mkdir -p %s'%out_dir)
	os.system('mkdir -p %s/data' %out_dir)
	os.system('mkdir -p %s/img' %out_dir)

	prior_correction = True 

	# initialize variables
	N_bins = 401
	P_dMfbyMf_dchifbychif_joint = np.ones((N_bins, N_bins))

	p = plt.figure(figsize=(5.2,5))
	ax1 = p.add_subplot(221)
	ax2 = p.add_subplot(224)
	ax3 = p.add_subplot(223)
	 
	color = ['r', 'orange', 'g', 'b']
	label = ['GW150914', 'GW170104']

	for (i, post_loc) in enumerate(post_loc_list):

		# read the posterior data for this event
		P_dMfbyMf_dchifbychif = np.loadtxt('%s/P_dMfbyMf_dchifbychif.dat.gz'%post_loc)
		dMfbyMf_vec = np.loadtxt('%s/dMfbyMf_vec.dat.gz'%post_loc)
		dchifbychif_vec = np.loadtxt('%s/dchifbychif_vec.dat.gz'%post_loc)

		conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
		s1_v1v2 = conf_v1v2.height_from_level(0.68)
		s2_v1v2 = conf_v1v2.height_from_level(0.95)
		ax3.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif), levels=(s1_v1v2, s2_v1v2), linewidths=(0.5,1.5), colors=color[i])

		# divide the posteior by the prior in (delta Mf/Mf, delta af/af) so that we use uninformed prior at the end 
		if prior_correction == True:
			prior_file = '../data/Prior_dMfbyMf_dchifbychif_Mfmin_1.0_Mfmax_500.0_chifmin-1.0_chifmax1.0_dMfbyMflim2.0_dchifbychiflim2.0'
			f = gzip.open(prior_file+".pklz",'rb')
			P_dMfbyMf_dchifbychif_pr_interp_obj = pickle.load(f)
			P_dMfbyMf_dchifbychif_pr = P_dMfbyMf_dchifbychif_pr_interp_obj(dMfbyMf_vec, dchifbychif_vec)
			P_dMfbyMf_dchifbychif /= P_dMfbyMf_dchifbychif_pr

		# removing nans and infinities 
		P_dMfbyMf_dchifbychif[np.isnan(P_dMfbyMf_dchifbychif)] = 0.
		P_dMfbyMf_dchifbychif[np.isinf(P_dMfbyMf_dchifbychif)] = 0.

		# replace all zeros in the posterior vector by a small number. 
		# This is to avoid multiplication by zeros while combining the posteiror 
		zidx = np.where(P_dMfbyMf_dchifbychif == 0.)
		P_dMfbyMf_dchifbychif[zidx] = 1e-16

		# Joint Probability distribution computation
		P_dMfbyMf_dchifbychif_joint = P_dMfbyMf_dchifbychif_joint*P_dMfbyMf_dchifbychif
		dx = np.mean(np.diff(dMfbyMf_vec))
		dy = np.mean(np.diff(dchifbychif_vec))
		P_dMfbyMf_dchifbychif_joint /= np.sum(P_dMfbyMf_dchifbychif_joint) * dx * dy
		P_dMfbyMf_dchifbychif /= np.sum(P_dMfbyMf_dchifbychif) * dx * dy

		# Marginalization to one-dimensional joint_posteriors
		P_dMfbyMf_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=0) * dy
		P_dchifbychif_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=1) * dx
		P_dMfbyMf = np.sum(P_dMfbyMf_dchifbychif, axis=0) * dy
		P_dchifbychif = np.sum(P_dMfbyMf_dchifbychif, axis=1) * dx

		# normalisation of individual posteriors
		P_dMfbyMf_joint /= np.sum(P_dMfbyMf_joint) * dx
		P_dchifbychif_joint /= np.sum(P_dchifbychif_joint) * dy
		P_dMfbyMf /= np.sum(P_dMfbyMf) * dx
		P_dchifbychif /= np.sum(P_dchifbychif) * dy

		# Calculation of confidence levels in the 2D posterior 
		conf_v1v2_joint = confidence(P_dMfbyMf_dchifbychif_joint)
		s1_v1v2_joint = conf_v1v2_joint.height_from_level(0.68)
		s2_v1v2_joint = conf_v1v2_joint.height_from_level(0.95)

		# calcualte the confidence levels and intervals in the marginalized 1d posteriors (joint) 
		s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_conf_intervals_in_1d(P_dMfbyMf_joint, dMfbyMf_vec)
		s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_conf_intervals_in_1d(P_dchifbychif_joint, dchifbychif_vec)

		# calcualte the confidence levels and intervals in the marginalized 1d posteriors (this event) 
		s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_conf_intervals_in_1d(P_dMfbyMf, dMfbyMf_vec)
		s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_conf_intervals_in_1d(P_dchifbychif, dchifbychif_vec)

		ax1.plot(dMfbyMf_vec, P_dMfbyMf, color=color[i], label=label[i]); ax2.plot(P_dchifbychif, dchifbychif_vec, color=color[i]);
		p.hold(True)

		# compute the credible level of the GR value of delta_Mf/Mf, delta_af/af) 
		conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
		gr_height = P_dMfbyMf_dchifbychif[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dchifbychif_vec))]
		gr_credib_level = conf_v1v2.level_from_height(gr_height)
		print 'Credible level of the GR value (%s): %.1f%% '%(label[i], 100.*gr_credib_level)

	# compute the credible level of the GR value from the combined posteriors 
	conf_v1v2 = confidence(P_dMfbyMf_dchifbychif_joint)
	gr_height = P_dMfbyMf_dchifbychif_joint[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dchifbychif_vec))]
	gr_credib_level = conf_v1v2.level_from_height(gr_height)
	print 'Credible level of the GR value (combined): %.1f%% '%(100.*gr_credib_level)

	np.savetxt(out_dir+'/data/dMfbyMf_vec.dat.gz', dMfbyMf_vec)
	np.savetxt(out_dir+'/data/dchifbychif_vec.dat.gz', dchifbychif_vec)
	np.savetxt(out_dir+'/data/P_dMfbyMf_dchifbychi_joint.dat.gz', P_dMfbyMf_dchifbychif_joint)
	np.savetxt(out_dir+'/data/P_dMfbyMf_joint.dat.gz', P_dMfbyMf_joint)
	np.savetxt(out_dir+'/data/P_dchifbychif_joint.dat.gz', P_dchifbychif_joint)

	ax1.plot(dMfbyMf_vec, P_dMfbyMf_joint, 'k',  label='Joint'); ax2.plot(P_dchifbychif_joint, dchifbychif_vec, 'k')
	ax3.contour(dMfbyMf_vec,dchifbychif_vec,gf(P_dMfbyMf_dchifbychif_joint), levels=(s1_v1v2_joint,s2_v1v2_joint), linewidths=(0.5,1.5), colors='k', label='Joint')
	ax3.plot(0, 0, 'k+', mew=0.5)
	#ax1.legend(loc='best', frameon=False)
	ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
	ax2.set_xlabel('$P(\Delta a_f/a_f)$')
	ax1.set_ylabel('$P(\Delta M_f/M_f)$')
	ax3.set_xlabel('$\Delta M_f/M_f$')
	ax3.set_ylabel('$\Delta a_f/a_f$')
	plt.tight_layout()
	plt.savefig(out_dir+'/img/dMfbyMf_dchifbychi_joint.pdf')

