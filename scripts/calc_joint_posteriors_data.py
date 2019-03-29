"""
Compute the joint posterior on Delta Mf/Mf and Delta chif/chif from multiple simulated BBH events

Abhirup Ghosh, P. Ajith, 2015-11-27
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy
import scipy.ndimage.filters as filter
import bayesian as ba
import time
import os
import os.path
import commands
import imrtestgr as tgr
import pickle, gzip
from scipy import interpolate
import optparse as op
from matplotlib import rc
import matplotlib
matplotlib.rc('text.latex', preamble = '\usepackage{txfonts}')
rc_params = {'backend': 'ps',
             'axes.labelsize': 10,
             'axes.titlesize': 10,
             'font.size': 12,
             'legend.fontsize': 8,
             'xtick.labelsize': 13,
             'ytick.labelsize': 13,
	     'font.family': 'Times New Roman'
             }

plt.rcParams.update(rc_params)

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
def calc_cred_intervals_in_1d(P, x):

	# find the value of P corresponding to 50% and 9% confidence heights
	conf = confidence(P)
	P_s1 = conf.height_from_level(0.5)
	P_s2 = conf.height_from_level(0.9)

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
	parser.add_option("-o","--out-dir", dest="out_dir", help="output directory")
	parser.add_option("-n","--N-bins", dest="N_bins", help="number of bins")
	parser.add_option("-l","--loc-file", dest="postloc_list",help="file with imrtgr results locations")
	(options, args) = parser.parse_args()
	postloc_list = options.postloc_list
        

	f=open(postloc_list,'r')
	lines=f.readlines()
	GW_events=[];GW_colors=[];GW_postloc_list=[]
	for x in lines:
    	   GW_events.append(x.split(' ')[0])
    	   GW_colors.append(x.split(' ')[1])
           GW_postloc_list.append((x.split(' ')[2]).split('\n')[0])  #
        f.close()
        print GW_events


	print GW_colors
	print GW_postloc_list

        print "These are the locations of imrtgr results of GW events:", GW_postloc_list

	out_dir = options.out_dir
	N_bins = int(options.N_bins)

	os.system('mkdir -p %s'%out_dir)
	os.system('mkdir -p %s/data' %out_dir)
	os.system('mkdir -p %s/img' %out_dir)

	# initialize variables
	P_dMfbyMf_dchifbychif_joint = np.ones((N_bins, N_bins))

	p = plt.figure(figsize=(5,5))
	ax1 = plt.subplot2grid((3,3), (0,0), colspan=2)
	ax2 = plt.subplot2grid((3,3), (1,2), rowspan=2)
	ax3 = plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)

	# load the pickle file containing the interpolant of the prior in delta_Mf/Mf, delta_af/a
	prior_file = '../data/Prior_normbymean_bbh_average_fits_precessing_m12min5.0_m12max170.00_a12zmin0.00_a12zmax0.99_epsilon2.5_sigma2.5_101bins_m1m2a1a2'
	f = gzip.open(prior_file+".pklz",'rb')
	P_dMfbyMf_dchifbychif_pr_interp_obj = pickle.load(f)

	for (i_event, post_loc) in enumerate(GW_postloc_list):

		# read the posterior data for this event
		P_dMfbyMf_dchifbychif = np.loadtxt('%s/P_dMfbyMf_dchifbychif.dat.gz'%post_loc)
		dMfbyMf_vec = np.loadtxt('%s/dMfbyMf_vec.dat.gz'%post_loc)
		dchifbychif_vec = np.loadtxt('%s/dchifbychif_vec.dat.gz'%post_loc)
		P_dMfbyMf_dchifbychif_interp_obj = scipy.interpolate.interp2d(dMfbyMf_vec, dchifbychif_vec, P_dMfbyMf_dchifbychif, fill_value=0., bounds_error=False)	
		dMfbyMf_vec = np.linspace(-2.0, 2.0, N_bins)
		P_dMfbyMf_dchifbychif = P_dMfbyMf_dchifbychif_interp_obj(dMfbyMf_vec, dchifbychif_vec)

		dx = np.mean(np.diff(dMfbyMf_vec))
                dy = np.mean(np.diff(dchifbychif_vec))
		P_dMfbyMf_dchifbychif /= np.sum(P_dMfbyMf_dchifbychif) * dx * dy

		# compute the credible regions of this posterior
		conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
		s1_v1v2 = conf_v1v2.height_from_level(0.5)
		s2_v1v2 = conf_v1v2.height_from_level(0.9)

		P_dMfbyMf = np.sum(P_dMfbyMf_dchifbychif, axis=0) * dy
                P_dchifbychif = np.sum(P_dMfbyMf_dchifbychif, axis=1) * dx

		P_dMfbyMf /= np.sum(P_dMfbyMf) * dx
                P_dchifbychif /= np.sum(P_dchifbychif) * dy

		# calcualte the credible levels and intervals in the marginalized 1d posteriors (this event)
                s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_cred_intervals_in_1d(P_dMfbyMf, dMfbyMf_vec)
                s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_cred_intervals_in_1d(P_dchifbychif, dchifbychif_vec)

		gr_height = P_dMfbyMf_dchifbychif[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dchifbychif_vec))]
                gr_credib_level = conf_v1v2.level_from_height(gr_height)
                print 'Credible level of the GR value (%s): %.1f%% '%(GW_events[i_event], 100.*gr_credib_level)
                print '90%% credible intervals on dMfbyMfmean (%s): %2.3f, %2.3f '%(GW_events[i_event], left2_v1, right2_v1)
                print '90%% credible intervals on dchifbychifmean (%s): %2.3f, %2.3f '%(GW_events[i_event], left2_v2, right2_v2)

                ax1.plot(dMfbyMf_vec, P_dMfbyMf, color=GW_colors[i_event], label=GW_events[i_event])
                #ax1.axvline(x = left2_v1, color=GW_colors[i_event], lw=0.5)
                #ax1.axvline(x = right2_v1, color=GW_colors[i_event], lw=0.5)
                ax2.plot(P_dchifbychif, dchifbychif_vec, color=GW_colors[i_event]);
                #ax2.axhline(y = left2_v2, color=GW_colors[i_event], lw=0.5)
                #ax2.axhline(y = right2_v2, color=GW_colors[i_event], lw=0.5)
		ax3.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif), levels=(s2_v1v2,), linewidths=(1.5,0.5), colors=GW_colors[i_event])

		# evaluate the prior interpolant
		P_dMfbyMf_dchifbychif_pr = P_dMfbyMf_dchifbychif_pr_interp_obj(dMfbyMf_vec, dchifbychif_vec)

		# removing nans
		P_dMfbyMf_dchifbychif[np.isnan(P_dMfbyMf_dchifbychif)] = 0.
		P_dMfbyMf_dchifbychif[np.isinf(P_dMfbyMf_dchifbychif)] = 0.

		# replace all zeros in the posterior vector by a small number.
		# This is to avoid multiplication by zeros while combining the posteiror
		zidx = np.where(P_dMfbyMf_dchifbychif == 0.)
		P_dMfbyMf_dchifbychif[zidx] = 1e-16

		if GW_events[i_event] == "LVT170729*":
			print "%s: not used for joint posterior computation"%GW_events[i_event]

		else:
			print "%s: used for joint posterior computation"%GW_events[i_event]
			print i_event
			# Joint Probability distribution computation. From the second event onwards, divide the posterior by the corresponding
			# prior so that we are multiplying the previous posterior with the likelihood of the current event
			if i_event == 0:
				P_dMfbyMf_dchifbychif_joint = P_dMfbyMf_dchifbychif_joint*P_dMfbyMf_dchifbychif
			elif i_event > 0:
				P_dMfbyMf_dchifbychif_joint = P_dMfbyMf_dchifbychif_joint*P_dMfbyMf_dchifbychif/P_dMfbyMf_dchifbychif_pr

			# removing nans and inf
			P_dMfbyMf_dchifbychif_joint[np.isnan(P_dMfbyMf_dchifbychif_joint)] = 0.
			P_dMfbyMf_dchifbychif_joint[np.isinf(P_dMfbyMf_dchifbychif_joint)] = 0.
			P_dMfbyMf_dchifbychif_joint /= np.sum(P_dMfbyMf_dchifbychif_joint) * dx * dy


	# Marginalization to one-dimensional prior posteriors
        P_dMfbyMf_pr = np.sum(P_dMfbyMf_dchifbychif_pr, axis=0) * dy
        P_dchifbychif_pr = np.sum(P_dMfbyMf_dchifbychif_pr, axis=1) * dx

        # normalisation of prior posteriors
        P_dMfbyMf_pr /= np.sum(P_dMfbyMf_pr) * dx
        P_dchifbychif_pr /= np.sum(P_dchifbychif_pr) * dy

	# Calculation of credible levels in the 2D posterior
        conf_v1v2_pr = confidence(P_dMfbyMf_dchifbychif_pr)
        s1_v1v2_pr = conf_v1v2_pr.height_from_level(0.5)
        s2_v1v2_pr = conf_v1v2_pr.height_from_level(0.9)

	# compute the credible level of the GR value from the combined posteriors
	conf_v1v2_joint = confidence(P_dMfbyMf_dchifbychif_joint)
	s1_v1v2_joint = conf_v1v2_joint.height_from_level(0.5)
        s2_v1v2_joint = conf_v1v2_joint.height_from_level(0.9)

	# Marginalization to one-dimensional joint_posteriors
        P_dMfbyMf_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=0) * dy
        P_dchifbychif_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=1) * dx

        # normalisation of individual posteriors
        P_dMfbyMf_joint /= np.sum(P_dMfbyMf_joint) * dx
        P_dchifbychif_joint /= np.sum(P_dchifbychif_joint) * dy

        # calcualte the credible levels and intervals in the marginalized 1d posteriors (joint)
        s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_cred_intervals_in_1d(P_dMfbyMf_joint, dMfbyMf_vec)
        s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_cred_intervals_in_1d(P_dchifbychif_joint, dchifbychif_vec)

	gr_height_joint = P_dMfbyMf_dchifbychif_joint[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dchifbychif_vec))]
	gr_credib_level_joint = conf_v1v2_joint.level_from_height(gr_height_joint)
	print 'Credible level of the GR value (combined): %.1f%% '%(100.*gr_credib_level_joint)
	print '90%% credible intervals on dMfbyMfmean (combined): %2.3f, %2.3f '%(left2_v1_joint, right2_v1_joint)
        print '90%% credible intervals on dchifbychifmean (combined): %2.3f, %2.3f '%(left2_v2_joint, right2_v2_joint)

	np.savetxt(out_dir+'/data/dMfbyMf_vec.dat.gz', dMfbyMf_vec)
	np.savetxt(out_dir+'/data/dchifbychif_vec.dat.gz', dchifbychif_vec)
	np.savetxt(out_dir+'/data/P_dMfbyMf_dchifbychi_joint.dat.gz', P_dMfbyMf_dchifbychif_joint)
	np.savetxt(out_dir+'/data/P_dMfbyMf_joint.dat.gz', P_dMfbyMf_joint)
	np.savetxt(out_dir+'/data/P_dchifbychif_joint.dat.gz', P_dchifbychif_joint)

	ax1.fill(dMfbyMf_vec, P_dMfbyMf_joint, 'darkred',  label='Joint', ls='--', alpha=0.5)
	ax1.plot(dMfbyMf_vec, P_dMfbyMf_pr, 'k',  lw=0.5, label='Prior', ls='--')
	#ax1.axvline(x = left2_v1_joint, color='k', lw=0.5)
        #ax1.axvline(x = right2_v1_joint, color='k', lw=0.5)
	ax2.fill(P_dchifbychif_joint, dchifbychif_vec, 'darkred', ls='--', alpha=0.5)
	ax2.plot(P_dchifbychif_pr, dchifbychif_vec, 'k', lw=0.5, ls='--')
	#ax2.axhline(y = left2_v2_joint, color='k', lw=0.5)
        #ax2.axhline(y = right2_v2_joint, color='k', lw=0.5)
	ax3.contourf(dMfbyMf_vec,dchifbychif_vec,gf(P_dMfbyMf_dchifbychif_joint), levels=(s2_v1v2_joint, np.inf), linewidths=(1.5,0.5), colors='darkred', label='Joint', linestyles=('--'), alpha=0.5)
	ax3.contour(dMfbyMf_vec,dchifbychif_vec,gf(P_dMfbyMf_dchifbychif_pr), levels=(s2_v1v2_pr,), linewidths=(0.5,0.5), colors='k', label='Prior', linestyles=('--'))
	ax3.plot(0, 0, 'k+', mew=1)
	ax1.set_xlim(-1.5,1.5)
	ax1.set_ylim(0,7)
	ax1.xaxis.tick_top()
	ax1.set_xticks(np.arange(-1.5, 1.5, 0.5))
        ax1.set_yticks(np.arange(0, 7, 1.))
	ax2.set_ylim(-1,1)
	ax2.set_xlim(0,7)
        ax2.yaxis.tick_right()
        ax2.set_xticks(np.arange(0, 8, 1.))
	ax2.set_yticks(np.arange(-1, 1.5, 0.5))
	ax3.set_xlim(-1.5,1.5)
        ax3.set_ylim(-1,1)
	ax3.set_xticks(np.arange(-1.5, 1.5, 0.5))
        ax3.set_yticks(np.arange(-1, 1.5, 0.5))
	ax1.legend(bbox_to_anchor=(1.05, 1.45), loc=2, borderaxespad=0.,frameon=False, fontsize="small")
	ax2.set_xlabel('$P(\Delta a_{\mathrm{f}}/\\bar{a}_{\mathrm{f}})$', fontsize='large', labelpad=1)
	ax1.set_ylabel('$P(\Delta M_{\mathrm{f}}/\\bar{M}_{\mathrm{f}})$', fontsize='large', labelpad=14)
	ax3.set_xlabel('$\Delta M_{\mathrm{f}}/\\bar{M}_{\mathrm{f}}$', fontsize='large', labelpad=1)
	ax3.set_ylabel('$\Delta a_{\mathrm{f}}/\\bar{a}_{\mathrm{f}}$', fontsize='large', labelpad=-1)
	plt.subplots_adjust(wspace=0., hspace=0.)
	#plt.tight_layout()
	plt.savefig(out_dir+'/img/dMfbyMf_dchifbychi_joint.pdf')
	plt.savefig(out_dir+'/img/dMfbyMf_dchifbychi_joint.png', dpi=300)
	plt.close()

