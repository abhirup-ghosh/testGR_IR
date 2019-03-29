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
postloc_lalinf_insp = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'
postloc_lalinf_ring = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/post-inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat'

postloc_tdimr_imr = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/IMR/emcee_samples.dat'
postloc_tdimr_insp = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/inspiral_132Hz_m0p004185/emcee_samples.dat'
postloc_tdimr_ring = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/post-inspiral_132Hz_m0p004185/emcee_samples.dat'

data = np.genfromtxt(postloc_lalinf_ring, dtype=None, names=True)
mc1, q1, dL1, i1, t01, phi01, ra1, sin_dec1, pol1 = data['mc'], data['q'], data['distance'], data['theta_jn'], data['time'], data['phase'], data['ra'], np.sin(data['dec']), data['psi']
m11, m21 = gw.comp_from_mcq(mc1, q1)
mf1, af1 = tgr.calc_final_mass_spin(m11, m21, np.zeros(len(m11)), np.zeros(len(m21)), 'nonprecspin_Healy2014')
Fp1,Fc1 = detector.overhead_antenna_pattern(ra1, np.arcsin(sin_dec1), pol1)

P_m11m21, m11_bins, m21_bins = np.histogram2d(m11, m21, bins=50, normed=True)
P_mc1q1, mc1_bins, q1_bins = np.histogram2d(mc1, q1, bins=50, normed=True)
P_mf1af1, mf1_bins, af1_bins = np.histogram2d(mf1, af1, bins=50, normed=True)

P_m11m21 = P_m11m21.T
P_mc1q1 = P_mc1q1.T
P_mf1af1 = P_mf1af1.T

conf_m11m21 = confidence(P_m11m21)
s1_m11m21 = conf_m11m21.height_from_level(0.5)
s2_m11m21 = conf_m11m21.height_from_level(0.9)

conf_mc1q1 = confidence(P_mc1q1)
s1_mc1q1 = conf_mc1q1.height_from_level(0.5)
s2_mc1q1 = conf_mc1q1.height_from_level(0.9)

conf_mf1af1 = confidence(P_mf1af1)
s1_mf1af1 = conf_mf1af1.height_from_level(0.5)
s2_mf1af1 = conf_mf1af1.height_from_level(0.9)

nwalkers, num_iter, ndim, n_burnin = 100, 20000, 9, 10000
mc2, q2, dL2, i2, t02, phi02, ra2, sin_dec2, pol2 = res.read_emcee_samples_9dim(postloc_tdimr_ring, nwalkers, num_iter, ndim, n_burnin)
m12, m22 = gw.comp_from_mcq(mc2, q2)
#mf2, af2 = tgr.calc_final_mass_spin(m12, m22, np.zeros(len(m12)), np.zeros(len(m22)), 'nonprecspin_Healy2014')
mf2, af2 = np.loadtxt('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/post-inspiral_132Hz_m0p004185/mfaf_samples.dat', unpack=True)
#np.savetxt('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/post-inspiral_132Hz_m0p004185/mfaf_samples.dat', np.c_[mf2, af2], header='mf af')
Fp2,Fc2 = detector.overhead_antenna_pattern(ra2, np.arcsin(sin_dec2), pol2)

P_m12m22, m12_bins, m22_bins = np.histogram2d(m12, m22, bins=50, normed=True)
P_mc2q2, mc2_bins, q2_bins = np.histogram2d(mc2, q2, bins=50, normed=True)
P_mf2af2, mf2_bins, af2_bins = np.histogram2d(mf2, af2, bins=50, normed=True)

P_m12m2 = P_m12m22.T
P_mc2q2 = P_mc2q2.T
P_mf2af2 = P_mf2af2.T

conf_m12m22 = confidence(P_m12m22)
s1_m12m22 = conf_m12m22.height_from_level(0.5)
s2_m12m22 = conf_m12m22.height_from_level(0.9)

conf_mc2q2 = confidence(P_mc2q2)
s1_mc2q2 = conf_mc2q2.height_from_level(0.5)
s2_mc2q2 = conf_mc2q2.height_from_level(0.9)
 
conf_mf2af2 = confidence(P_mf2af2)
s1_mf2af2 = conf_mf2af2.height_from_level(0.5)
s2_mf2af2 = conf_mf2af2.height_from_level(0.9)

fig = plt.figure(figsize=(6,3))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.contour(mc1_bins[:-1], q1_bins[:-1], tgr.gf(P_mc1q1), levels=(s2_mc1q1,s1_mc1q1), linewidths=(1,1.5), linestyles=('solid', 'dashed'),colors='red')
ax1.contour(mc2_bins[:-1], q2_bins[:-1], tgr.gf(P_mc2q2), levels=(s2_mc2q2,s1_mc2q2), linewidths=(1,1.5), linestyles=('solid', 'dashed'),colors='k')
ax1.set_xlabel('$\mathcal{M}_c$')
ax1.set_ylabel('$q$')
ax2.contour(mf1_bins[:-1], af1_bins[:-1], tgr.gf(P_mf1af1), levels=(s2_mf1af1,s1_mf1af1), linewidths=(1,1.5), linestyles=('solid', 'dashed'),colors='red')
ax2.contour(mf2_bins[:-1], af2_bins[:-1], tgr.gf(P_mf2af2), levels=(s2_mf2af2,s1_mf2af2), linewidths=(1,1.5), linestyles=('solid', 'dashed'),colors='k')
ax2.set_xlabel('$M_f$')
ax2.set_ylabel('$a_f$')
ax2.set_xlim([min(min(mf1), min(mf2)),70])
ax2.set_ylim([0.5,max(max(af1), max(af2))])
plt.tight_layout()
plt.savefig('/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_reruns_3629/2D_comparison_ring.png')
