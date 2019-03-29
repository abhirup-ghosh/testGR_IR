import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import imrtgrutils_final as tgr
import read_emcee_samples as res
import standard_gwtransf as gw
import scipy
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

post_loc_LAL_IMR = "/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/IMR/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat"
post_loc_LAL_insp = "/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat"
post_loc_LAL_ring = "/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalinference/post-inspiral/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat"
post_loc_TD_IMR = "/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_50000iter/emcee_samples.dat"
post_loc_TD_insp = "/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_50000iter_inspiral/emcee_samples.dat" 
post_loc_TD_ring = "/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/PE_exercises/9_param_lalsim_TD_losc_t0inwhiten_100000iter_post-inspiral_m0p00927/emcee_samples.dat"

fit_formula = "bbh_average_fits_precessing"


#data_LAL_IMR = np.genfromtxt(post_loc_LAL_IMR, dtype=None, names=True)
#m1_LAL_IMR, m2_LAL_IMR = data_LAL_IMR['m1'], data_LAL_IMR['m2']
#chi1_LAL_IMR, chi2_LAL_IMR, chi1z_LAL_IMR, chi2z_LAL_IMR, tilt1_LAL_IMR, tilt2_LAL_IMR, phi12_LAL_IMR = np.zeros(len(m1_LAL_IMR)), np.zeros(len(m1_LAL_IMR)), np.zeros(len(m1_LAL_IMR)), np.zeros(len(m1_LAL_IMR)), np.zeros(len(m1_LAL_IMR)), np.zeros(len(m1_LAL_IMR)), np.zeros(len(m1_LAL_IMR))
Mf_LAL_IMR, chif_LAL_IMR = np.loadtxt('tdimr_old_scripts/mfchif_samples_LAL_IMR.dat', unpack=True)#tgr.calc_final_mass_spin(m1_LAL_IMR, m2_LAL_IMR, chi1_LAL_IMR, chi2_LAL_IMR, chi1z_LAL_IMR, chi2z_LAL_IMR, tilt1_LAL_IMR, tilt2_LAL_IMR, phi12_LAL_IMR, fit_formula)


#data_LAL_insp = np.genfromtxt(post_loc_LAL_insp, dtype=None, names=True)
#m1_LAL_insp, m2_LAL_insp = data_LAL_insp['m1'], data_LAL_insp['m2']
#chi1_LAL_insp, chi2_LAL_insp, chi1z_LAL_insp, chi2z_LAL_insp, tilt1_LAL_insp, tilt2_LAL_insp, phi12_LAL_insp = np.zeros(len(m1_LAL_insp)), np.zeros(len(m1_LAL_insp)), np.zeros(len(m1_LAL_insp)), np.zeros(len(m1_LAL_insp)), np.zeros(len(m1_LAL_insp)), np.zeros(len(m1_LAL_insp)), np.zeros(len(m1_LAL_insp))
Mf_LAL_insp, chif_LAL_insp =  np.loadtxt('tdimr_old_scripts/mfchif_samples_LAL_insp.dat', unpack=True)#tgr.calc_final_mass_spin(m1_LAL_insp, m2_LAL_insp, chi1_LAL_insp, chi2_LAL_insp, chi1z_LAL_insp, chi2z_LAL_insp, tilt1_LAL_insp, tilt2_LAL_insp, phi12_LAL_insp, fit_formula)

Mf_LAL_ring, chif_LAL_ring =  np.loadtxt('tdimr_old_scripts/mfchif_samples_LAL_ring.dat', unpack=True)

#nwalkers, num_iter, ndim, n_burnin = 100, 50000, 9, 30000
#mc_TD_IMR, q_TD_IMR, dL_TD_IMR, i_TD_IMR, t0_TD_IMR, phi0_TD_IMR, ra_TD_IMR, sin_dec_TD_IMR, pol_TD_IMR = res.read_emcee_samples_9dim(post_loc_TD_IMR, nwalkers, num_iter, ndim, n_burnin)
#m1_TD_IMR, m2_TD_IMR = gw.comp_from_mcq(mc_TD_IMR, q_TD_IMR)
Mf_TD_IMR, chif_TD_IMR = np.loadtxt('tdimr_old_scripts/mfchif_samples_TD_IMR_burnin30000.dat', unpack=True)#tgr.calc_final_mass_spin(m1_TD_IMR, m2_TD_IMR, np.zeros(len(m1_TD_IMR)), np.zeros(len(m1_TD_IMR)), np.zeros(len(m1_TD_IMR)), np.zeros(len(m1_TD_IMR)), np.zeros(len(m1_TD_IMR)), np.zeros(len(m1_TD_IMR)), np.zeros(len(m1_TD_IMR)), fit_formula)

#nwalkers, num_iter, ndim, n_burnin = 100, 50000, 9, 30000
#mc_TD_insp, q_TD_insp, dL_TD_insp, i_TD_insp, t0_TD_insp, phi0_TD_insp, ra_TD_insp, sin_dec_TD_insp, pol_TD_insp = res.read_emcee_samples_9dim(post_loc_TD_insp, nwalkers, num_iter, ndim, n_burnin)
#m1_TD_insp, m2_TD_insp = gw.comp_from_mcq(mc_TD_insp, q_TD_insp)
Mf_TD_insp, chif_TD_insp = np.loadtxt('tdimr_old_scripts/mfchif_samples_TD_insp_burnin30000.dat', unpack=True)#tgr.calc_final_mass_spin(m1_TD_insp, m2_TD_insp, np.zeros(len(m1_TD_insp)), np.zeros(len(m1_TD_insp)), np.zeros(len(m1_TD_insp)), np.zeros(len(m1_TD_insp)), np.zeros(len(m1_TD_insp)), np.zeros(len(m1_TD_insp)), np.zeros(len(m1_TD_insp)), fit_formula)

Mf_TD_ring, chif_TD_ring = np.loadtxt('tdimr_old_scripts/mfchif_samples_TD_ring_burnin60000.dat', unpack=True)#tgr.calc_final_mass_spin(m1_TD_ring, m2_TD_ring, np.zeros(len(m1_TD_ring)), np.zeros(len(m1_TD_ring)), np.zeros(len(m1_TD_ring)), np.zeros(len(m1_TD_ring)), np.zeros(len(m1_TD_ring)), np.zeros(len(m1_TD_ring)), np.zeros(len(m1_TD_ring)), fit_formula)

print '... started plotting'

N_bins = 1601

Mf_max = max(max(Mf_LAL_IMR), max(Mf_TD_IMR), max(Mf_LAL_insp), max(Mf_TD_insp), max(Mf_LAL_ring), max(Mf_TD_ring))
Mf_min = min(min(Mf_LAL_IMR), min(Mf_TD_IMR), min(Mf_LAL_insp), min(Mf_TD_insp), min(Mf_LAL_ring), min(Mf_TD_ring))
chif_max = max(max(chif_LAL_IMR), max(chif_TD_IMR), max(chif_LAL_insp), max(chif_TD_insp), max(chif_LAL_ring), max(chif_TD_ring))
chif_min = min(min(chif_LAL_IMR), min(chif_TD_IMR), min(chif_LAL_insp), min(chif_TD_insp), min(chif_LAL_ring), min(chif_TD_ring))

Mf_bins = np.linspace(Mf_min, Mf_max, N_bins)
chif_bins = np.linspace(chif_min, chif_max, N_bins)

dMf = np.mean(np.diff(Mf_bins))
dchif = np.mean(np.diff(chif_bins))

Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
chif_intp = (chif_bins[:-1] + chif_bins[1:])/2.

P_Mfchif_LAL_IMR, Mf_bins, chif_bins = np.histogram2d(Mf_LAL_IMR, chif_LAL_IMR, bins=(Mf_bins, chif_bins), normed=True)
P_Mfchif_LAL_insp, Mf_bins, chif_bins = np.histogram2d(Mf_LAL_insp, chif_LAL_insp, bins=(Mf_bins, chif_bins), normed=True)
P_Mfchif_LAL_ring, Mf_bins, chif_bins = np.histogram2d(Mf_LAL_ring, chif_LAL_ring, bins=(Mf_bins, chif_bins), normed=True)

P_Mfchif_TD_IMR, Mf_bins, chif_bins = np.histogram2d(Mf_TD_IMR, chif_TD_IMR, bins=(Mf_bins, chif_bins), normed=True)
P_Mfchif_TD_insp, Mf_bins, chif_bins = np.histogram2d(Mf_TD_insp, chif_TD_insp, bins=(Mf_bins, chif_bins), normed=True)
P_Mfchif_TD_ring, Mf_bins, chif_bins = np.histogram2d(Mf_TD_ring, chif_TD_ring, bins=(Mf_bins, chif_bins), normed=True)
  

P_Mfchif_LAL_IMR = P_Mfchif_LAL_IMR.T
P_Mfchif_LAL_insp = P_Mfchif_LAL_insp.T
P_Mfchif_LAL_ring = P_Mfchif_LAL_ring.T

P_Mfchif_TD_IMR = P_Mfchif_TD_IMR.T
P_Mfchif_TD_insp = P_Mfchif_TD_insp.T
P_Mfchif_TD_ring = P_Mfchif_TD_ring.T

conf_Mfchif_LAL_IMR = confidence(P_Mfchif_LAL_IMR)
s1_Mfchif_LAL_IMR = conf_Mfchif_LAL_IMR.height_from_level(0.5)
s2_Mfchif_LAL_IMR = conf_Mfchif_LAL_IMR.height_from_level(0.9)

conf_Mfchif_LAL_insp = confidence(P_Mfchif_LAL_insp)
s1_Mfchif_LAL_insp = conf_Mfchif_LAL_insp.height_from_level(0.5)
s2_Mfchif_LAL_insp = conf_Mfchif_LAL_insp.height_from_level(0.9)

conf_Mfchif_LAL_ring = confidence(P_Mfchif_LAL_ring)
s1_Mfchif_LAL_ring = conf_Mfchif_LAL_ring.height_from_level(0.5)
s2_Mfchif_LAL_ring = conf_Mfchif_LAL_ring.height_from_level(0.9)

conf_Mfchif_TD_IMR = confidence(P_Mfchif_TD_IMR)
s1_Mfchif_TD_IMR = conf_Mfchif_TD_IMR.height_from_level(0.5)
s2_Mfchif_TD_IMR = conf_Mfchif_TD_IMR.height_from_level(0.9)

conf_Mfchif_TD_insp = confidence(P_Mfchif_TD_insp)
s1_Mfchif_TD_insp = conf_Mfchif_TD_insp.height_from_level(0.5)
s2_Mfchif_TD_insp = conf_Mfchif_TD_insp.height_from_level(0.9)

conf_Mfchif_TD_ring = confidence(P_Mfchif_TD_ring)
s1_Mfchif_TD_ring = conf_Mfchif_TD_ring.height_from_level(0.5)
s2_Mfchif_TD_ring = conf_Mfchif_TD_ring.height_from_level(0.9)

plt.figure(figsize=(5,5))
plt.subplot2grid((3,3), (0,0), colspan=2)
#plt.hist(Mf_LAL_IMR, bins=51, histtype='step', color='k', ls='dashed',normed=True)
#plt.hist(Mf_LAL_insp, bins=51, histtype='step', color='orange', ls='dashed',normed=True)
plt.hist(Mf_LAL_ring, bins=101, histtype='step', color='r', ls='dashed',normed=True)
#plt.hist(Mf_TD_IMR, bins=51, histtype='step', color='k',normed=True)
#plt.hist(Mf_TD_insp, bins=51, histtype='step', color='orange',normed=True)
plt.hist(Mf_TD_ring, bins=101, histtype='step', color='k',normed=True)
#plt.xlabel('$M_f$')
#plt.xlim([60,64])
plt.xticks([])
plt.ylabel('$P(M_f)$')
#plt.grid()

plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
#CS_LAL_IMR = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_LAL_IMR), levels=(s2_Mfchif_LAL_IMR,), linewidths=(1,), linestyles=('dashed',),colors='k')
#CS_LAL_insp = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_LAL_insp), levels=(s2_Mfchif_LAL_insp,), linewidths=(1,), linestyles=('dashed',),colors='orange')
CS_LAL_ring = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_LAL_ring), levels=(s2_Mfchif_LAL_ring,), linewidths=(1,), linestyles=('dashed',),colors='red')
#CS_TD_IMR = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_TD_IMR), levels=(s2_Mfchif_TD_IMR,), linewidths=(1,), colors='k')
#CS_TD_insp = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_TD_insp), levels=(s2_Mfchif_TD_insp,), linewidths=(1,), colors='orange')
CS_TD_ring = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_TD_ring), levels=(s2_Mfchif_TD_ring,), linewidths=(1,), colors='k')
#plt.axvline(x=Mf_inj, ls='--', color='k')
#plt.axhline(y=chif_inj, ls='--', color='k')
#plt.xlim([min(np.append(Mf_LAL_IMR, Mf_TD_IMR)), max(np.append(Mf_LAL_IMR, Mf_TD_IMR))])
#plt.ylim([min(np.append(chif_LAL_IMR, chif_TD_IMR)), max(np.append(chif_LAL_IMR, chif_TD_IMR))])
plt.xlabel('$M_f~[M_\odot]$')
plt.ylabel('$\chi_f$')
plt.grid()
#plt.xlim([60,64])
#plt.ylim([0.64,chif_max])
plt.grid()

plt.subplot2grid((3,3), (1,2), rowspan=2)
#plt.hist(chif_LAL_IMR, bins=51, normed=True, histtype='step', color='k', ls='dashed', orientation='horizontal')
#plt.hist(chif_LAL_insp, bins=51, normed=True, histtype='step', color='orange',ls='dashed',orientation='horizontal')
plt.hist(chif_LAL_ring, bins=101, normed=True, histtype='step', color='r',ls='dashed',orientation='horizontal')
#plt.hist(chif_TD_IMR, bins=51, normed=True, histtype='step', color='k',orientation='horizontal')
#plt.hist(chif_TD_insp, bins=51, normed=True, histtype='step', color='orange',orientation='horizontal')
plt.hist(chif_TD_ring, bins=101, normed=True, histtype='step', color='k',orientation='horizontal')
#plt.ylabel('$\chi_f$')
plt.ylim([0.64,chif_max])
plt.yticks([])
plt.xlabel('$P(\chi_f)$')
#plt.grid()
plt.tight_layout()
plt.savefig('./TDFD_comparison_ring.png')
