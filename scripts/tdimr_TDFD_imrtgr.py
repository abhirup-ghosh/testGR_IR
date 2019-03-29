import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import read_emcee_samples as res
import pycbc
from pycbc  import  detector
import standard_gwtransf as gw
import imrtestgr as tgr
import scipy
import scipy.signal as ss
from scipy import interpolate

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

postloc_lalinf_imr = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/IMR/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'
postloc_lalinf_i = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'
postloc_lalinf_r = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/post-inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'

postloc_tdimr_imr = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/IMR/emcee_samples.dat'
postloc_tdimr_i = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/inspiral_132Hz_m0p004185/emcee_samples.dat'
postloc_tdimr_r = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/post-inspiral_132Hz_m0p004185/emcee_samples.dat'

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)

data_i = np.genfromtxt(postloc_lalinf_i, dtype=None, names=True)
mc_i, q_i, dL_i, i_i, t0_i, phi0_i, ra_i, sin_dec_i, pol_i = data_i['mc'], data_i['q'], data_i['distance'], data_i['theta_jn'], data_i['time'], data_i['phase'], data_i['ra'], np.sin(data_i['dec']), data_i['psi']
m1_i, m2_i = gw.comp_from_mcq(mc_i, q_i)
mf_i, af_i = tgr.calc_final_mass_spin(m1_i, m2_i, np.zeros(len(m1_i)), np.zeros(len(m2_i)), 'nonprecspin_Healy2014')
P_mfaf_i, mf_i_bins, af_i_bins = np.histogram2d(mf_i, af_i, bins=50, normed=True)
P_mfaf_i = P_mfaf_i.T
conf_mfaf_i = confidence(P_mfaf_i)
s1_mfaf_i = conf_mfaf_i.height_from_level(0.5)
s2_mfaf_i = conf_mfaf_i.height_from_level(0.9)

data_r = np.genfromtxt(postloc_lalinf_r, dtype=None, names=True)
mc_r, q_r, dL_r, i_r, t0_r, phi0_r, ra_r, sin_dec_r, pol_r = data_r['mc'], data_r['q'], data_r['distance'], data_r['theta_jn'], data_r['time'], data_r['phase'], data_r['ra'], np.sin(data_r['dec']), data_r['psi']
m1_r, m2_r = gw.comp_from_mcq(mc_r, q_r)
mf_r, af_r = tgr.calc_final_mass_spin(m1_r, m2_r, np.zeros(len(m1_r)), np.zeros(len(m2_r)), 'nonprecspin_Healy2014')
P_mfaf_r, mf_r_bins, af_r_bins = np.histogram2d(mf_r, af_r, bins=50, normed=True)
P_mfaf_r = P_mfaf_r.T
conf_mfaf_r = confidence(P_mfaf_r)
s1_mfaf_r = conf_mfaf_r.height_from_level(0.5)
s2_mfaf_r = conf_mfaf_r.height_from_level(0.9)

data_imr = np.genfromtxt(postloc_lalinf_imr, dtype=None, names=True)
mc_imr, q_imr, dL_imr, i_imr, t0_imr, phi0_imr, ra_imr, sin_dec_imr, pol_imr = data_imr['mc'], data_imr['q'], data_imr['distance'], data_imr['theta_jn'], data_imr['time'], data_imr['phase'], data_imr['ra'], np.sin(data_imr['dec']), data_imr['psi']
m1_imr, m2_imr = gw.comp_from_mcq(mc_imr, q_imr)
mf_imr, af_imr = tgr.calc_final_mass_spin(m1_imr, m2_imr, np.zeros(len(m1_imr)), np.zeros(len(m2_imr)), 'nonprecspin_Healy2014')
P_mfaf_imr, mf_imr_bins, af_imr_bins = np.histogram2d(mf_imr, af_imr, bins=50, normed=True)
P_mfaf_imr = P_mfaf_imr.T
conf_mfaf_imr = confidence(P_mfaf_imr)
s1_mfaf_imr = conf_mfaf_imr.height_from_level(0.5)
s2_mfaf_imr = conf_mfaf_imr.height_from_level(0.9)


ax.contour(mf_i_bins[:-1], af_i_bins[:-1], tgr.gf(P_mfaf_i), levels=(s2_mfaf_i,), linewidths=(1,1.5), linestyles=('solid', 'solid'),colors='orange')
ax.contour(mf_r_bins[:-1], af_r_bins[:-1], tgr.gf(P_mfaf_r), levels=(s2_mfaf_r,), linewidths=(1,1.5), linestyles=('solid', 'solid'),colors='red')
ax.contour(mf_imr_bins[:-1], af_imr_bins[:-1], tgr.gf(P_mfaf_imr), levels=(s2_mfaf_imr,), linewidths=(1,1.5), linestyles=('solid', 'solid'),colors='k')


nwalkers, num_iter, ndim, n_burnin = 100, 20000, 9, 10000
mc_i, q_i, dL_i, i_i, t0_i, phi0_i, ra_i, sin_dec_i, pol_i = res.read_emcee_samples_9dim(postloc_tdimr_i, nwalkers, num_iter, ndim, n_burnin)
m1_i, m2_i = gw.comp_from_mcq(mc_i, q_i)
mf_i, af_i = np.loadtxt('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/inspiral_132Hz_m0p004185/mfaf_samples.dat', unpack=True)
P_mfaf_i, mf_i_bins, af_i_bins = np.histogram2d(mf_i, af_i, bins=50, normed=True)
P_mfaf_i = P_mfaf_i.T
conf_mfaf_i = confidence(P_mfaf_i)
s1_mfaf_i = conf_mfaf_i.height_from_level(0.5)
s2_mfaf_i = conf_mfaf_i.height_from_level(0.9)

nwalkers, num_iter, ndim, n_burnin = 100, 20000, 9, 10000
mc_r, q_r, dL_r, i_r, t0_r, phi0_r, ra_r, sin_dec_r, pol_r = res.read_emcee_samples_9dim(postloc_tdimr_r, nwalkers, num_iter, ndim, n_burnin)
m1_r, m2_r = gw.comp_from_mcq(mc_r, q_r)
mf_r, af_r = np.loadtxt('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/post-inspiral_132Hz_m0p004185/mfaf_samples.dat', unpack=True)
P_mfaf_r, mf_r_bins, af_r_bins = np.histogram2d(mf_r, af_r, bins=50, normed=True)
P_mfaf_r = P_mfaf_r.T
conf_mfaf_r = confidence(P_mfaf_r)
s1_mfaf_r = conf_mfaf_r.height_from_level(0.5)
s2_mfaf_r = conf_mfaf_r.height_from_level(0.9)

nwalkers, num_iter, ndim, n_burnin = 100, 20000, 9, 10000
mc_imr, q_imr, dL_imr, i_imr, t0_imr, phi0_imr, ra_imr, sin_dec_imr, pol_imr = res.read_emcee_samples_9dim(postloc_tdimr_imr, nwalkers, num_iter, ndim, n_burnin)
m1_imr, m2_imr = gw.comp_from_mcq(mc_imr, q_imr)
mf_imr, af_imr = np.loadtxt('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/IMR/mfaf_samples.dat', unpack=True)
P_mfaf_imr, mf_imr_bins, af_imr_bins = np.histogram2d(mf_imr, af_imr, bins=50, normed=True)
P_mfaf_imr = P_mfaf_imr.T
conf_mfaf_imr = confidence(P_mfaf_imr)
s1_mfaf_imr = conf_mfaf_imr.height_from_level(0.5)
s2_mfaf_imr = conf_mfaf_imr.height_from_level(0.9)

ax.contour(mf_i_bins[:-1], af_i_bins[:-1], tgr.gf(P_mfaf_i), levels=(s2_mfaf_i,), linewidths=(1,1.5), linestyles=('dashed', 'dashed'),colors='orange')
ax.contour(mf_r_bins[:-1], af_r_bins[:-1], tgr.gf(P_mfaf_r), levels=(s2_mfaf_r,), linewidths=(1,1.5), linestyles=('dashed', 'dashed'),colors='red')
ax.contour(mf_imr_bins[:-1], af_imr_bins[:-1], tgr.gf(P_mfaf_imr), levels=(s2_mfaf_imr,), linewidths=(1,1.5), linestyles=('dashed', 'dashed'),colors='k')
ax.set_xlim([55,70])
ax.set_ylim([0.55,0.68])
ax.set_xlabel('$M_f$')
ax.set_ylabel('$a_f$')
plt.tight_layout()
plt.savefig('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/imrtgr_comparison.png')
