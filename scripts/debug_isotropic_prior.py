#import matplotlib as mpl
#mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import lalinference.imrtgr.imrtgrutils as tgr
import scipy
import scipy.signal as ss
from scipy import interpolate
from optparse import OptionParser
import pickle, gzip

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

# location of IMRPPv2 posterior samples
post_loc = '/home/abhirup/Documents/Work/O2/2017/January/04/1167559936p5991/G268556/lalinference/20170217_prodruns_officialbuild_d5de7def/imrtgr_results/IMRPhenomPv2_nonprecspin_Healy2014_2017-02-21_fhigh_insp143Hz_flow_ring_143Hz_401bins_pc_debug/data'
#post_loc_seob = '/home/abhirup.ghosh/Documents/Work/O2/2017/January/04/1167559936p5991/G268556/lalinference/20170217_prodruns_officialbuild_d5de7def/imrtgr_results/SEOBNRv4_ROM_nonprecspin_Healy2014_2017-02-21_fhigh_insp143Hz_flow_ring_143Hz_401bins_pc_debug/data'

# loading the data from the posterior location
Mf, chif = np.loadtxt('%s/Mfchif.dat.gz'%post_loc)
Mf = (Mf[:-1] + Mf[1:])/2.
chif = (chif[:-1] + chif[1:])/2.
dMf = np.mean(np.diff(Mf))
dchif = np.mean(np.diff(chif))

P_Mfchif_i_pre = np.loadtxt('%s/P_Mfchif_i_pre.dat.gz'%post_loc) # imr consistency test output before prior correction
P_Mfchif_i_pre /= np.sum(P_Mfchif_i_pre) * dMf * dchif
conf_Mfchif_i_pre = confidence(P_Mfchif_i_pre)
s1_Mfchif_i_pre = conf_Mfchif_i_pre.height_from_level(0.68)
s2_Mfchif_i_pre = conf_Mfchif_i_pre.height_from_level(0.95)

P_Mfchif_i = np.loadtxt('%s/P_Mfchif_i.dat.gz'%post_loc) # imr consistency test output after prior correction
conf_Mfchif_i = confidence(P_Mfchif_i)
s1_Mfchif_i = conf_Mfchif_i.height_from_level(0.68)
s2_Mfchif_i = conf_Mfchif_i.height_from_level(0.95)


# loading the different priors
prior_lalinference = '/home/abhirup/opt/lalsuite_lalinference_o2_d5de7def/share/lalinference/imrtgr_prior_data/Prior_Mfaf_nonprec_Healy2014_M1-500_isotropic.pklz'
prior_recomp = '/home/abhirup/Documents/Work/testGR_IR/data/Prior_Mfaf_nonprecspin_Healy2014_comp_mass_min1.0_comp_mass_max500.0_comp_spin_min0.0_comp_spin_max1.0_isotropic.pklz'
prior_short = '/home/abhirup/Documents/Work/testGR_IR/data/Prior_Mfaf_nonprecspin_Healy2014_comp_mass_min5.0_comp_mass_max170.0_comp_spin_min0.0_comp_spin_max1.0_isotropic.pklz'

f1 = gzip.open(prior_lalinference,'rb')
P_Mfchif_pr_lalinf_interp_obj = pickle.load(f1)
P_Mfchif_pr_lalinf = P_Mfchif_pr_lalinf_interp_obj(Mf, chif)

f2 = gzip.open(prior_recomp,'rb')
P_Mfchif_pr_recomp_interp_obj = pickle.load(f2)
P_Mfchif_pr_recomp = P_Mfchif_pr_recomp_interp_obj(Mf, chif)

f3 = gzip.open(prior_short,'rb')
P_Mfchif_pr_short_interp_obj = pickle.load(f3)
P_Mfchif_pr_short = P_Mfchif_pr_short_interp_obj(Mf, chif)

# prior correction using lalinference prior
P_Mfchif_i_post_lalinf = P_Mfchif_i_pre / P_Mfchif_pr_lalinf
P_Mfchif_i_post_lalinf[np.isinf(P_Mfchif_i_post_lalinf)] = 0.
P_Mfchif_i_post_lalinf[np.isnan(P_Mfchif_i_post_lalinf)] = 0.
conf_Mfchif_i_post_lalinf = confidence(P_Mfchif_i_post_lalinf)
s1_Mfchif_i_post_lalinf = conf_Mfchif_i_post_lalinf.height_from_level(0.68)
s2_Mfchif_i_post_lalinf = conf_Mfchif_i_post_lalinf.height_from_level(0.95)

# prior correction using recomputed prior
P_Mfchif_i_post_recomp = P_Mfchif_i_pre / P_Mfchif_pr_recomp
P_Mfchif_i_post_recomp[np.isinf(P_Mfchif_i_post_recomp)] = 0.
P_Mfchif_i_post_recomp[np.isnan(P_Mfchif_i_post_recomp)] = 0.
conf_Mfchif_i_post_recomp = confidence(P_Mfchif_i_post_recomp)
s1_Mfchif_i_post_recomp = conf_Mfchif_i_post_recomp.height_from_level(0.68)
s2_Mfchif_i_post_recomp = conf_Mfchif_i_post_recomp.height_from_level(0.95)

# prior correction using short prior
P_Mfchif_i_post_short = P_Mfchif_i_pre / P_Mfchif_pr_short
P_Mfchif_i_post_short[np.isinf(P_Mfchif_i_post_short)] = 0.
P_Mfchif_i_post_short[np.isnan(P_Mfchif_i_post_short)] = 0.
conf_Mfchif_i_post_short = confidence(P_Mfchif_i_post_short)
s1_Mfchif_i_post_short = conf_Mfchif_i_post_short.height_from_level(0.68)
s2_Mfchif_i_post_short = conf_Mfchif_i_post_short.height_from_level(0.95)


plt.figure(figsize=(16,16))
plt.subplot(331)
plt.pcolormesh(Mf, chif, P_Mfchif_i_pre, cmap='YlOrBr')
plt.colorbar()
plt.contour(Mf, chif, tgr.gf(P_Mfchif_i_pre), levels=(s1_Mfchif_i_pre,), linewidths=(1,1.5))
plt.xlim([0, 100])
plt.ylim([0,1])
plt.title('pre prior corr')
plt.subplot(332)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_lalinf, cmap='YlOrBr')
plt.colorbar()
plt.xlim([0, 100])
plt.ylim([0,1])
plt.title('$P_{lalinference}(Mf,chif)$')
plt.subplot(333)
plt.pcolormesh(Mf, chif, P_Mfchif_i_post_lalinf, cmap='YlOrBr')
plt.colorbar()
plt.contour(Mf, chif, tgr.gf(P_Mfchif_i_post_lalinf), levels=(s1_Mfchif_i_post_lalinf,), linewidths=(1,1.5))
plt.xlim([0, 100])
plt.ylim([0,1])
plt.title('post prior corr')

plt.subplot(335)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_recomp, cmap='YlOrBr')
plt.colorbar()
plt.xlim([0, 100])
plt.ylim([0,1])
plt.title('$P_{recomp}(Mf,chif)$')
plt.subplot(336)
plt.pcolormesh(Mf, chif, P_Mfchif_i_post_recomp, cmap='YlOrBr')
plt.colorbar()
plt.contour(Mf, chif, tgr.gf(P_Mfchif_i_post_recomp), levels=(s1_Mfchif_i_post_recomp,), linewidths=(1,1.5))
plt.xlim([0, 100])
plt.ylim([0,1])
plt.title('post prior corr')

plt.subplot(338)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_short, cmap='YlOrBr')
plt.colorbar()
plt.xlim([0, 100])
plt.ylim([0,1])
plt.title('$P_{short}(Mf,chif)$')
plt.subplot(339)
plt.pcolormesh(Mf, chif, P_Mfchif_i_post_short, cmap='YlOrBr')
plt.colorbar()
plt.contour(Mf, chif, tgr.gf(P_Mfchif_i_post_short), levels=(s1_Mfchif_i_post_short,), linewidths=(1,1.5))
plt.xlim([0, 100])
plt.ylim([0,1])
plt.tight_layout()
#plt.savefig('debug_imrpp.png')

plt.figure(figsize=(16,8))
plt.subplot(231)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_lalinf, cmap='YlOrBr')
plt.colorbar()
plt.xlim([0, 100])
plt.ylim([0,1])
plt.title('$P_{lalinference}(Mf,chif)$')
plt.subplot(232)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_recomp, cmap='YlOrBr')
plt.colorbar()
plt.title('$P_{recomp}(Mf,chif)$')
plt.xlim([0, 100])
plt.ylim([0,1])
plt.subplot(233)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_short, cmap='YlOrBr')
plt.colorbar()
plt.title('$P_{short}(Mf,chif)$')
plt.xlim([0, 100])
plt.ylim([0,1])
plt.subplot(234)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_lalinf - P_Mfchif_pr_recomp, cmap='YlOrBr')
plt.title('$P_{lalinference}(Mf,chif) - P_{recomp}(Mf,chif)$')
plt.colorbar()
plt.subplot(235)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_lalinf - P_Mfchif_pr_short, cmap='YlOrBr')
plt.colorbar()
plt.title('$P_{lalinference}(Mf,chif) - P_{short}(Mf,chif)$')
plt.subplot(236)
plt.pcolormesh(Mf, chif, P_Mfchif_pr_recomp - P_Mfchif_pr_short, cmap='YlOrBr')
plt.colorbar()
plt.title('$P_{recomp}(Mf,chif) - P_{short}(Mf,chif)$')
plt.show()
