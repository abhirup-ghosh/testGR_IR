import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import scipy.ndimage.filters as filter
import bayesian as ba

#import common
#from common import events, event_colors
#from common import sample_files
#from common import contour_dir
#from common import rc_params

# Make the font match the document
rc_params = {'backend': 'ps',
             'axes.labelsize': 10,
             'axes.titlesize': 10,
             'font.size': 12,
             'legend.fontsize': 10,
             'xtick.labelsize': 13,
             'ytick.labelsize': 13,
             'font.family': 'Times New Roman'#,
             #'font.family': 'sans-serif',
             #'font.sans-serif': ['Bitstream Vera Sans']
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

	post_loc = '/Users/nkjm/testing_GR/testGR_IR/data_from_Alice/GW170104/imrtgr_results/IMRPhenomPv2_bbh_average_fits_precessing_2017-03-31_fhigh_insp143Hz_flow_ring_143Hz_401bins_evolved_samples_normbymean_final_noMfafprior/data'

	P_dMfbyMf_dchifbychif = np.loadtxt('%s/P_dMfbyMf_dchifbychif.dat.gz'%post_loc)
	P_dMfbyMf  = np.loadtxt('%s/P_dMfbyMf.dat.gz'%post_loc)
	P_dchifbychif  = np.loadtxt('%s/P_dchifbychif.dat.gz'%post_loc)
	dMfbyMf_vec = np.loadtxt('%s/dMfbyMf_vec.dat.gz'%post_loc)
	dchifbychif_vec = np.loadtxt('%s/dchifbychif_vec.dat.gz'%post_loc)

	p = plt.figure(figsize=(5,5))
	ax1 = plt.subplot2grid((3,3), (0,0), colspan=2)
	ax2 = plt.subplot2grid((3,3), (1,2), rowspan=2)
	ax3 = plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)

	conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
	s1_v1v2 = conf_v1v2.height_from_level(0.5)
	s2_v1v2 = conf_v1v2.height_from_level(0.9)
	ax3.pcolormesh(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif), cmap='Purples')
	ax3.contour(dMfbyMf_vec,dchifbychif_vec, gf(P_dMfbyMf_dchifbychif), levels=(s2_v1v2, s1_v1v2), linewidths=(1.5,0.5), colors='k')


	s1_v1, s2_v1, left1_v1, right1_v1, left2_v1, right2_v1 = calc_cred_intervals_in_1d(P_dMfbyMf, dMfbyMf_vec)
	s1_v2, s2_v2, left1_v2, right1_v2, left2_v2, right2_v2 = calc_cred_intervals_in_1d(P_dchifbychif, dchifbychif_vec)

	ax1.plot(dMfbyMf_vec, P_dMfbyMf, color='k')
	ax1.fill_between(dMfbyMf_vec, P_dMfbyMf, 0, color='k', alpha=0.2)
	ax1.axvline(x = left2_v1, color='k', lw=0.5, ls='--')
	ax1.axvline(x = right2_v1, color='k', lw=0.5, ls='--')
	ax2.plot(P_dchifbychif, dchifbychif_vec, color='k');
	ax2.fill_betweenx(dchifbychif_vec, P_dchifbychif, 0, color='k', alpha=0.2)
	ax2.axhline(y = left2_v2, color='k', lw=0.5, ls='--')
	ax2.axhline(y = right2_v2, color='k', lw=0.5, ls='--')

	ax1.set_xlim(-1,1)
	ax1.set_ylim(0,3)
	ax2.set_ylim(-1,1)
	ax2.set_xlim(0,3)
	ax2.set_xlabel('$P(\Delta a_{\mathrm{f}}/\\bar{a}_{\mathrm{f}})$', fontsize='large', labelpad=1)
	ax1.set_ylabel('$P(\Delta M_{\mathrm{f}}/\\bar{M}_{\mathrm{f}})$', fontsize='large', labelpad=14)
	ax3.set_xlabel('$\Delta M_{\mathrm{f}}/\\bar{M}_{\mathrm{f}}$', fontsize='large', labelpad=1)
	ax3.set_ylabel('$\Delta a_{\mathrm{f}}/\\bar{a}_{\mathrm{f}}$', fontsize='large', labelpad=-1)
	ax1.set_xticks(np.arange(-1, 1.5, 0.5))
	ax1.xaxis.tick_top()
	ax1.set_yticks(np.arange(0, 4, 1.))
	ax2.set_yticks(np.arange(-1, 1.5, 0.5))
	ax2.yaxis.tick_right()
	ax2.set_xticks(np.arange(0, 4, 1.))
	ax3.set_xticks(np.arange(-1, 1., 0.5))
	ax3.set_yticks(np.arange(-1, 1.5, 0.5))
	plt.subplots_adjust(wspace=0., hspace=0.)

	plt.savefig('dMfbyMf_dchifbychi_GW170104.pdf')
        plt.savefig('dMfbyMf_dchifbychi_GW170104.png', dpi=300)
        plt.close()
