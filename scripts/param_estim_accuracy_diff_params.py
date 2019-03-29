import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import lalinference.imrtgr.imrtgrutils as tgr
import numpy as np
import scipy
import scipy.signal as ss
from scipy import interpolate
from optparse import OptionParser
import time, os
import pickle, gzip
import sys
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
from pylal import git_version
import scipy.ndimage.filters as filter
import standard_gwtransf as gw
#import plotsettings 

from matplotlib import rc
import matplotlib

def gf(P):
	return filter.gaussian_filter(P, sigma=2.0)

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

injection_no = 39

#GR locations
data_i = np.genfromtxt('/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-12-25/uniform_compmass_spins/injection_%d/inspiral/1126285216-0/H1L1/posterior_samples.dat'%injection_no,dtype=None, names=True)#np.genfromtxt('/home/archis/public_html/testGR_IR/runs/simulations_modGR/2015-11-27/modGR_INJECTED-1126285214-4_q1_GR_new_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_aLIGO_EARLY_HIGH_seed19076/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/inspiral/1126285216.0-0/H1L1/posterior_samples.dat',dtype=None, names=True)

data_r = np.genfromtxt('/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-12-25/uniform_compmass_spins/injection_%d/post-inspiral/1126285216-0/H1L1/posterior_samples.dat'%injection_no,dtype=None, names=True)#np.genfromtxt('/home/archis/public_html/testGR_IR/runs/simulations_modGR/2015-11-27/modGR_INJECTED-1126285214-4_q1_GR_new_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_aLIGO_EARLY_HIGH_seed19076/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/post-inspiral/1126285216.0-0/H1L1/posterior_samples.dat',dtype=None, names=True)

data_imr = np.genfromtxt('/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-12-25/uniform_compmass_spins/injection_%d/IMR/1126285216-0/H1L1/posterior_samples.dat'%injection_no,dtype=None, names=True)#np.genfromtxt('/home/archis/public_html/testGR_IR/runs/simulations_modGR/2015-11-27/modGR_INJECTED-1126285214-4_q1_GR_new_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_aLIGO_EARLY_HIGH_seed19076/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/IMR/1126285216.0-0/H1L1/posterior_samples.dat',dtype=None, names=True)

#modGR locations
#data_i = np.genfromtxt('/home/archis/public_html/testGR_IR/runs/simulations_modGR/2015-12-07/modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_aLIGO_EARLY_HIGH_seed22347/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/inspiral/1126285216.0-0/H1L1/posterior_samples.dat',dtype=None, names=True)

#data_r = np.genfromtxt('/home/archis/public_html/testGR_IR/runs/simulations_modGR/2015-12-07/modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_aLIGO_EARLY_HIGH_seed22347/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/post-inspiral/1126285216.0-0/H1L1/posterior_samples.dat',dtype=None, names=True)

#data_imr = np.genfromtxt('/home/archis/public_html/testGR_IR/runs/simulations_modGR/2015-12-07/modGR_INJECTED-1126285214-4_q1_a2_400_r0_32_M70_dist420_incl2p54_ra1p76_decm1p23_psi1p6_flow20_aLIGO_EARLY_HIGH_seed22347/SEOBNRv2_ROM_DoubleSpinthreePointFivePN_seglen4/IMR/1126285216.0-0/H1L1/posterior_samples.dat',dtype=None, names=True)

#injected values
fit_formula = 'nonprecspin_Healy2014'

m1_inj, m2_inj, chi1_inj, chi2_inj = np.genfromtxt('/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2015-12-25/uniform_compmass_spins/injection_%d/IMR/injection.txt'%injection_no, usecols=(1,0,6,7))

print m1_inj, m2_inj, chi1_inj, chi2_inj

mtotal_inj = m1_inj + m2_inj
q_inj = m2_inj/m1_inj
eta_inj = gw.eta_from_q(q_inj)
chi_eff_inj = 0.
mc_inj = gw.mc_from_comp(m1_inj, m2_inj)
Mf_inj, chif_inj = tgr.calc_final_mass_spin(m1_inj, m2_inj, chi1_inj, chi2_inj, fit_formula)
out_dir = '/home/abhirup/public_html/imrtestgr/ER8/G184098/investigations/posteriors_on_diff_params/2016-02-06'
os.system('mkdir -p %s'%out_dir)

x1_list = ['mtotal', 'mf', 'mc']
x2_list = ['chi_eff', 'af', 'q']
x1_inj_list = [mtotal_inj, Mf_inj, mc_inj]
x2_inj_list = [1., 1., 1.]
x1_label_list = ['$M/M_{inj}$', '$M_f/M_{f, inj}$', '$M_c/M_{c, inj}$']
x2_label_list = ['$\chi _{eff}$', '$\chi_f$', '$q$']
x2_lim_list = [-1., 0., 0.]

idx = 1

plt.figure(figsize=(15,15))
plt.text(.9, 0.7, 'inspiral', color='orange')
plt.text(.9, 0.6, 'post-inspiral', color='red')
plt.text(.9, 0.5, 'IMR', color='k')
for (x2, x2_inj, x2_label, x2_lim) in zip(x2_list, x2_inj_list, x2_label_list, x2_lim_list):
  for (x1, x1_inj, x1_label) in zip(x1_list, x1_inj_list, x1_label_list):
    print idx
    x1_i, x2_i = data_i[x1]/x1_inj, data_i[x2]/x2_inj
    x1_r, x2_r = data_r[x1]/x1_inj, data_r[x2]/x2_inj
    x1_imr, x2_imr = data_imr[x1]/x1_inj, data_imr[x2]/x2_inj

    P_i, x1_i_bins, x2_i_bins = np.histogram2d(x1_i, x2_i, bins=51, normed=True); P_i = P_i.T; conf_i = confidence(P_i); s1_i = conf_i.height_from_level(0.68); s2_i = conf_i.height_from_level(0.95)
    P_r, x1_r_bins, x2_r_bins = np.histogram2d(x1_r, x2_r, bins=51, normed=True); P_r = P_r.T; conf_r = confidence(P_r); s1_r = conf_r.height_from_level(0.68); s2_r = conf_r.height_from_level(0.95)
    P_imr, x1_imr_bins, x2_imr_bins = np.histogram2d(x1_imr, x2_imr, bins=51, normed=True); P_imr = P_imr.T; conf_imr = confidence(P_imr); s1_imr = conf_imr.height_from_level(0.68); s2_imr = conf_imr.height_from_level(0.95)

    plt.subplot(3,3,idx)
    plt.contour(x1_i_bins[:-1], x2_i_bins[:-1], gf(P_i), levels=(s1_i,), colors='orange', linewidths=(1,1.5), label='inspiral')
    plt.contour(x1_r_bins[:-1], x2_r_bins[:-1], gf(P_r), levels=(s1_r,), colors='r', linewidths=(1,1.5))
    plt.contour(x1_imr_bins[:-1], x2_imr_bins[:-1], gf(P_imr), levels=(s1_imr,), colors='k', linewidths=(1,1.5))
    plt.grid()
    plt.xlim(0.7, 1.2)
    plt.ylim(x2_lim,1.)
    plt.xlabel('%s'%x1_label)
    plt.ylabel('%s'%x2_label)
    plt.hold(True)    
    idx = idx + 1
plt.tight_layout()
plt.savefig('%s/IMR_overlap_plot_alt_param_GRinj_spinning_inj_%d.png'%(out_dir, injection_no), dpi=300)
#plt.savefig('%s/IMR_overlap_plot_alt_param_modGRinj_a2_400.png'%out_dir, dpi=300)
