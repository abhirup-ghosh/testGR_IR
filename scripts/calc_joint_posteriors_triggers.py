import sys, lal, os, commands, numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
import lalinference.imrtgr.imrtgrutils as tgr
import bayesian as ba

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

def set_tick_sizes(ax, major, minor):
  for l in ax.get_xticklines() + ax.get_yticklines():
    l.set_markersize(major)
  for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
    tick.tick1line.set_markersize(minor)
    tick.tick2line.set_markersize(minor)
  ax.xaxis.LABELPAD=10.
  ax.xaxis.OFFSETTEXTPAD=10.

#Module for confidence calculations
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

#IMR_consistency_test output locations
gw1_loc = '/home/abhirup/public_html/imrtestgr/ER8/G184098/SEOBNRv2_IMRPhenomPv2_combined_nonprecspin_Healy2014_2016-04-20_fhigh_insp132Hz_flow_ring_132Hz_prodruns_C01_o1_lalinference_20160402'
gw2_loc = '/home/abhirup/public_html/imrtestgr/O1/G211117/IMRPhenomPv2_nonprecspin_Healy2014_2016-04-20_fhigh_insp100Hz_flow_ring_100Hz_prod_runs_C01_SpinFix_wo_priorcorr_401bins'


P_dMfbyMf_dchifbychif_1 = np.loadtxt('%s/data/P_dMfbyMf_dchifbychif.dat.gz'%gw1_loc); dMfbyMf_vec_1 = np.loadtxt('%s/data/dMfbyMf_vec.dat.gz'%gw1_loc); dchifbychif_vec_1 = np.loadtxt('%s/data/dchifbychif_vec.dat.gz'%gw1_loc)

P_dMfbyMf_dchifbychif_2 = np.loadtxt('%s/data/P_dMfbyMf_dchifbychif.dat.gz'%gw2_loc); dMfbyMf_vec_2 = np.loadtxt('%s/data/dMfbyMf_vec.dat.gz'%gw2_loc); dchifbychif_vec_2 = np.loadtxt('%s/data/dchifbychif_vec.dat.gz'%gw2_loc)

if np.array_equal(dMfbyMf_vec_1, dMfbyMf_vec_2) == True and np.array_equal(dchifbychif_vec_1, dchifbychif_vec_2) == True:
	dMfbyMf_vec = dMfbyMf_vec_1
	dchifbychif_vec = dchifbychif_vec_1

# Joint Probability distribution computation
P_dMfbyMf_dchifbychif_joint = P_dMfbyMf_dchifbychif_1 * P_dMfbyMf_dchifbychif_2
dx = np.mean(np.diff(dMfbyMf_vec))
dy = np.mean(np.diff(dchifbychif_vec))
P_dMfbyMf_dchifbychif_joint /= np.sum(P_dMfbyMf_dchifbychif_joint) * dx * dy    # normalization

# Marginalization to one-dimensional joint_posteriors
P_dMfbyMf_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=0) * dy
P_dchifbychif_joint = np.sum(P_dMfbyMf_dchifbychif_joint, axis=1) * dx

# Calculation of confidence levels in the 2D posterior 
s1_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint, 0.68)
s2_v1v2_joint = ba.nsigma_value(P_dMfbyMf_dchifbychif_joint, 0.95)

# calcualte the confidence levels and intervals in the marginalized 1d posteriors (joint) 
s1_v1_joint, s2_v1_joint, left1_v1_joint, right1_v1_joint, left2_v1_joint, right2_v1_joint = calc_conf_intervals_in_1d(P_dMfbyMf_joint, dMfbyMf_vec)
s1_v2_joint, s2_v2_joint, left1_v2_joint, right1_v2_joint, left2_v2_joint, right2_v2_joint = calc_conf_intervals_in_1d(P_dchifbychif_joint, dchifbychif_vec)

#plotting
plt.figure(figsize=(5,5))
plt.subplot2grid((3,3), (0,0), colspan=2)
plt.plot(dMfbyMf_vec, P_dMfbyMf_joint, color='k', lw=1)
plt.axvline(x=left1_v1_joint, color='c', lw=0.5, ls='-.')
plt.axvline(x=right1_v1_joint, color='c', lw=0.5, ls='-.')
plt.axvline(x=left2_v1_joint, color='b', lw=0.5, ls='-.')
plt.axvline(x=right2_v1_joint, color='b', lw=0.5, ls='-.')
#plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$P(\Delta M_f/M_f)$')
#plt.grid()

plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
plt.pcolormesh(dMfbyMf_vec,dchifbychif_vec,P_dMfbyMf_dchifbychif_joint, cmap='YlOrBr')
plt.contour(dMfbyMf_vec,dchifbychif_vec,tgr.gf(P_dMfbyMf_dchifbychif_joint), levels=(s1_v1v2_joint,s2_v1v2_joint), linewidths=(1,1.5))
plt.plot(0, 0, 'k+', ms=12, mew=2)
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta \chi_f/\chi_f$')
plt.xlim([-1.,1.])
plt.ylim([-1.,1.])
plt.grid()

plt.subplot2grid((3,3), (1,2), rowspan=2)
plt.plot(P_dchifbychif_joint, dchifbychif_vec,'k', lw=1)
plt.axhline(y=left1_v2_joint, color='c', lw=0.5, ls='-.')
plt.axhline(y=right1_v2_joint, color='c', lw=0.5, ls='-.')
plt.axhline(y=left2_v2_joint, color='b', lw=0.5, ls='-.')
plt.axhline(y=right2_v2_joint, color='b', lw=0.5, ls='-.')
#plt.ylabel('$\Delta \chi_f/\chi_f$')
plt.xlabel('$P(\Delta \chi_f/\chi_f)$')
#plt.grid()
plt.savefig('dMfbyMfdchifbychif.png', dpi=300)
