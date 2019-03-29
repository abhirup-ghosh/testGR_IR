import sys, lal, os, commands, numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python'))
import lalinference.imrtgr.imrtgrutils as tgr

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

n_bins = 401

spin_fit_formula = 'nonprecspin_Healy2014'
test_no = 'G184098'#'G268556' 
insp_fhigh = 132
ring_flow = 132
root_tag = 'prodruns_C01_o2_lalinference_20170224'
waveform = 'SEOBNRv4_IMRPhenomPv2_combined'
date = '2017-02-24'

tag = '%s_%s_%s_fhigh_insp%dHz_flow_ring_%dHz_%s' %(waveform, spin_fit_formula, date, insp_fhigh, ring_flow, root_tag)
out_dir = '/home/abhirup.ghosh/Documents/Work/ER8/2015/September/15/1126259462p3910/G184098/lalinference/20170217_prodruns_officialbuild_d5de7def/imrtgr_results/%s' %tag
os.system('mkdir -p %s'%out_dir)
os.system('mkdir -p %s/data' %out_dir)
os.system('mkdir -p %s/img' %out_dir)

#IMR_consistency_test output locations
post_loc_imrpp = '/home/abhirup.ghosh/Documents/Work/ER8/2015/September/15/1126259462p3910/G184098/lalinference/20170217_prodruns_officialbuild_d5de7def/imrtgr_results/IMRPhenomPv2_nonprecspin_Healy2014_2017-02-23_fhigh_insp132Hz_flow_ring_132Hz_401bins_pc_final/data'
post_loc_seob = '/home/abhirup.ghosh/Documents/Work/ER8/2015/September/15/1126259462p3910/G184098/lalinference/20170217_prodruns_officialbuild_d5de7def/imrtgr_results/SEOBNRv4_ROM_nonprecspin_Healy2014_2017-02-23_fhigh_insp132Hz_flow_ring_132Hz_401bins_pc_final/data'

#IMRPhenomPv2 posteriors
Mf_imrpp, chif_imrpp = np.loadtxt('%s/Mfchif.dat.gz'%post_loc_imrpp)
Mf_imrpp = (Mf_imrpp[:-1] + Mf_imrpp[1:])/2.
chif_imrpp = (chif_imrpp[:-1] + chif_imrpp[1:])/2.
dMf_imrpp = np.mean(np.diff(Mf_imrpp))
dchif_imrpp = np.mean(np.diff(chif_imrpp))

P_Mfchif_i_imrpp = np.loadtxt('%s/P_Mfchif_i.dat.gz'%post_loc_imrpp)
P_Mfchif_r_imrpp = np.loadtxt('%s/P_Mfchif_r.dat.gz'%post_loc_imrpp)
P_Mfchif_imr_imrpp = np.loadtxt('%s/P_Mfchif_imr.dat.gz'%post_loc_imrpp)

dMfbyMf_vec_imrpp = np.loadtxt('%s/dMfbyMf_vec.dat.gz'%post_loc_imrpp)
dchifbychif_vec_imrpp = np.loadtxt('%s/dchifbychif_vec.dat.gz'%post_loc_imrpp)
P_dMfbyMf_dchifbychif_imrpp = np.loadtxt('%s/P_dMfbyMf_dchifbychif.dat.gz'%post_loc_imrpp)

#SEOBNRv2 posteriors
Mf_seob, chif_seob = np.loadtxt('%s/Mfchif.dat.gz'%post_loc_seob)
Mf_seob = (Mf_seob[:-1] + Mf_seob[1:])/2.
chif_seob = (chif_seob[:-1] + chif_seob[1:])/2.
dMf_seob = np.mean(np.diff(Mf_seob))
dchif_seob = np.mean(np.diff(chif_seob))

P_Mfchif_i_seob = np.loadtxt('%s/P_Mfchif_i.dat.gz'%post_loc_seob)
P_Mfchif_r_seob = np.loadtxt('%s/P_Mfchif_r.dat.gz'%post_loc_seob)
P_Mfchif_imr_seob = np.loadtxt('%s/P_Mfchif_imr.dat.gz'%post_loc_seob)

dMfbyMf_vec_seob = np.loadtxt('%s/dMfbyMf_vec.dat.gz'%post_loc_seob)
dchifbychif_vec_seob = np.loadtxt('%s/dchifbychif_vec.dat.gz'%post_loc_seob)
P_dMfbyMf_dchifbychif_seob = np.loadtxt('%s/P_dMfbyMf_dchifbychif.dat.gz'%post_loc_seob)

#normalisation of inspiral, ringdown, and IMR posteriors
P_Mfchif_i_imrpp /= np.sum(P_Mfchif_i_imrpp) * dMf_imrpp * dchif_imrpp
P_Mfchif_r_imrpp /= np.sum(P_Mfchif_r_imrpp) * dMf_imrpp * dchif_imrpp
P_Mfchif_imr_imrpp /= np.sum(P_Mfchif_imr_imrpp) * dMf_imrpp * dchif_imrpp

P_Mfchif_i_seob /= np.sum(P_Mfchif_i_seob) * dMf_seob * dchif_seob
P_Mfchif_r_seob /= np.sum(P_Mfchif_r_seob) * dMf_seob * dchif_seob
P_Mfchif_imr_seob /= np.sum(P_Mfchif_imr_seob) * dMf_seob * dchif_seob

#check that (dMfbyMf_vec, dchifbychif_vec) is identically defined between (-1,1) for IMRPP and SEOB, so that any one of them can be selected as the x- and y-bins for plotting P(dMfbyMf_vec, dchifbychif_vec)
if np.array_equal(dMfbyMf_vec_imrpp, dMfbyMf_vec_seob) == True and np.array_equal(dchifbychif_vec_imrpp, dchifbychif_vec_seob) == True:
  dMfbyMf_vec = dMfbyMf_vec_imrpp
  dchifbychif_vec = dchifbychif_vec_imrpp

diff_dMfbyMf = np.mean(np.diff(dMfbyMf_vec))
diff_dchifbychif = np.mean(np.diff(dchifbychif_vec))

#-------------------------------------------------------------------------------
#defining a uniform grid for plotting IMR overlap plot
#-------------------------------------------------------------------------------

#defining Mf, chif limits
Mf_lim = max(max(Mf_imrpp), max(Mf_seob))
chif_lim = max(max(chif_imrpp), max(chif_seob))

#defining Mf, chif bins
Mf_bins = np.linspace(-Mf_lim, Mf_lim, n_bins)
chif_bins = np.linspace(-chif_lim, chif_lim, n_bins)

Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
chif_intp = (chif_bins[:-1] + chif_bins[1:])/2.
dMf_intp = np.mean(np.diff(Mf_intp))
dchif_intp = np.mean(np.diff(chif_intp))

#interpolating posteriors on a uniform grid defined by (Mf_intp, chif_intp)
P_Mfchif_i_imrpp_inter_object = scipy.interpolate.interp2d(Mf_imrpp, chif_imrpp, P_Mfchif_i_imrpp, fill_value=0., bounds_error=False)
P_Mfchif_r_imrpp_inter_object = scipy.interpolate.interp2d(Mf_imrpp, chif_imrpp, P_Mfchif_r_imrpp, fill_value=0., bounds_error=False)
P_Mfchif_imr_imrpp_inter_object = scipy.interpolate.interp2d(Mf_imrpp, chif_imrpp, P_Mfchif_imr_imrpp, fill_value=0., bounds_error=False)

P_Mfchif_i_seob_inter_object = scipy.interpolate.interp2d(Mf_seob, chif_seob, P_Mfchif_i_seob, fill_value=0., bounds_error=False)
P_Mfchif_r_seob_inter_object = scipy.interpolate.interp2d(Mf_seob, chif_seob, P_Mfchif_r_seob, fill_value=0., bounds_error=False)
P_Mfchif_imr_seob_inter_object = scipy.interpolate.interp2d(Mf_seob, chif_seob, P_Mfchif_imr_seob, fill_value=0., bounds_error=False)

P_Mfchif_i_imrpp = P_Mfchif_i_imrpp_inter_object(Mf_intp, chif_intp)
P_Mfchif_r_imrpp = P_Mfchif_r_imrpp_inter_object(Mf_intp, chif_intp)
P_Mfchif_imr_imrpp = P_Mfchif_imr_imrpp_inter_object(Mf_intp, chif_intp)

P_Mfchif_i_seob = P_Mfchif_i_seob_inter_object(Mf_intp, chif_intp)
P_Mfchif_r_seob = P_Mfchif_r_seob_inter_object(Mf_intp, chif_intp)
P_Mfchif_imr_seob = P_Mfchif_imr_seob_inter_object(Mf_intp, chif_intp)

#adding posteriors
P_Mfchif_i = P_Mfchif_i_imrpp + P_Mfchif_i_seob
P_Mfchif_r = P_Mfchif_r_imrpp + P_Mfchif_r_seob
P_Mfchif_imr = P_Mfchif_imr_imrpp + P_Mfchif_imr_seob
P_dMfbyMf_dchifbychif = P_dMfbyMf_dchifbychif_imrpp + P_dMfbyMf_dchifbychif_seob

# Marginalization to one-dimensional joint_posteriors
P_dMfbyMf = np.sum(P_dMfbyMf_dchifbychif, axis=0) * diff_dchifbychif
P_dchifbychif = np.sum(P_dMfbyMf_dchifbychif, axis=1) * diff_dMfbyMf

#normalisation of combined posteriors
P_Mfchif_i /= np.sum(P_Mfchif_i) * dMf_intp * dchif_intp
P_Mfchif_r /= np.sum(P_Mfchif_r) * dMf_intp * dchif_intp
P_Mfchif_imr /= np.sum(P_Mfchif_imr) * dMf_intp * dchif_intp
P_dMfbyMf_dchifbychif /= np.sum(P_dMfbyMf_dchifbychif) * np.mean(np.diff(dMfbyMf_vec)) * np.mean(np.diff(dchifbychif_vec))
P_dMfbyMf /= np.sum(P_dMfbyMf) * diff_dMfbyMf 
P_dchifbychif /= np.sum(P_dchifbychif) * diff_dchifbychif 

#calculating confidence levels
conf_Mfchif_i = confidence(P_Mfchif_i)
s1_Mfchif_i = conf_Mfchif_i.height_from_level(0.68)
s2_Mfchif_i = conf_Mfchif_i.height_from_level(0.95)

conf_Mfchif_r = confidence(P_Mfchif_r)
s1_Mfchif_r = conf_Mfchif_r.height_from_level(0.68)
s2_Mfchif_r = conf_Mfchif_r.height_from_level(0.95)

conf_Mfchif_imr = confidence(P_Mfchif_imr)
s1_Mfchif_imr = conf_Mfchif_imr.height_from_level(0.68)
s2_Mfchif_imr = conf_Mfchif_imr.height_from_level(0.95)

conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
s1_v1v2 = conf_v1v2.height_from_level(0.68)
s2_v1v2 = conf_v1v2.height_from_level(0.95)

conf_v1 = confidence(P_dMfbyMf)
s1_v1 = conf_v1.height_from_level(0.68)
s2_v1 = conf_v1.height_from_level(0.95)

conf_v2 = confidence(P_dchifbychif)
s1_v2 = conf_v2.height_from_level(0.68)
s2_v2 = conf_v2.height_from_level(0.95)

# Calculation of condifence edges
left1_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])
right1_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s1_v1)[0]])

left2_v1 = min(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])
right2_v1 = max(dMfbyMf_vec[np.where(P_dMfbyMf>=s2_v1)[0]])

left1_v2 = min(dchifbychif_vec[np.where(P_dchifbychif>s1_v2)[0]])
right1_v2 = max(dchifbychif_vec[np.where(P_dchifbychif>s1_v2)[0]])

left2_v2 = min(dchifbychif_vec[np.where(P_dchifbychif>s2_v2)[0]])
right2_v2 = max(dchifbychif_vec[np.where(P_dchifbychif>s2_v2)[0]])

#plotting
plt.figure(figsize=(5,5))
CSi = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_i), levels=(s1_Mfchif_i,s2_Mfchif_i), linewidths=(1,1.5), colors='orange')
CSr = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_r), levels=(s1_Mfchif_r,s2_Mfchif_r), linewidths=(1,1.5), colors='red')
CSimr = plt.contour(Mf_bins[:-1], chif_bins[:-1], tgr.gf(P_Mfchif_imr), levels=(s1_Mfchif_imr,s2_Mfchif_imr), linewidths=(1,1.5), colors='k')
plt.xlim([40,160])
plt.ylim([0,1])
plt.xlabel('$M_f~[M_\odot]$')
plt.ylabel('$\chi_f$')
plt.grid()
plt.savefig('%s/img/IMR_overlap.png'%out_dir, dpi=300)

plt.figure(figsize=(5,5))
plt.subplot2grid((3,3), (0,0), colspan=2)
plt.plot(dMfbyMf_vec, P_dMfbyMf, color='k', lw=1)
plt.axvline(x=left1_v1, color='c', lw=0.5, ls='-.')
plt.axvline(x=right1_v1, color='c', lw=0.5, ls='-.')
plt.axvline(x=left2_v1, color='b', lw=0.5, ls='-.')
plt.axvline(x=right2_v1, color='b', lw=0.5, ls='-.')
#plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$P(\Delta M_f/M_f)$')
#plt.grid()

plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
plt.pcolormesh(dMfbyMf_vec,dchifbychif_vec,P_dMfbyMf_dchifbychif, cmap='YlOrBr')
plt.contour(dMfbyMf_vec,dchifbychif_vec,tgr.gf(P_dMfbyMf_dchifbychif), levels=(s1_v1v2,s2_v1v2), linewidths=(1,1.5))
plt.plot(0, 0, 'k+', ms=12, mew=2)
plt.xlabel('$\Delta M_f/M_f$')
plt.ylabel('$\Delta \chi_f/\chi_f$')
plt.xlim([-1.,1.])
plt.ylim([-1.,1.])
plt.grid()

plt.subplot2grid((3,3), (1,2), rowspan=2)
plt.plot(P_dchifbychif, dchifbychif_vec,'k', lw=1)
plt.axhline(y=left1_v2, color='c', lw=0.5, ls='-.')
plt.axhline(y=right1_v2, color='c', lw=0.5, ls='-.')
plt.axhline(y=left2_v2, color='b', lw=0.5, ls='-.')
plt.axhline(y=right2_v2, color='b', lw=0.5, ls='-.')
#plt.ylabel('$\Delta \chi_f/\chi_f$')
plt.xlabel('$P(\Delta \chi_f/\chi_f)$')
#plt.grid()
plt.savefig('%s/img/dMfbyMfdchifbychif.png'%out_dir, dpi=300)

# compute the confidence region corresponding to the GR value (delta_Mf/Mf = 0, delta_chif/chif = 0). 
# the 'confidence' class is defined on top of this script 
conf_v1v2 = confidence(P_dMfbyMf_dchifbychif)
gr_height = P_dMfbyMf_dchifbychif[np.argmin(abs(dMfbyMf_vec)), np.argmin(abs(dchifbychif_vec))] # taking value closest to (0,0)
gr_conf_level = conf_v1v2.level_from_height(gr_height)
print '... no deviation from GR above %.1f%% confidence level'%(100.*gr_conf_level)


np.savetxt(out_dir+'/data/Mfchif.dat.gz', (Mf_bins,chif_bins))
np.savetxt(out_dir+'/data/P_Mfchif_i.dat.gz', P_Mfchif_i)
np.savetxt(out_dir+'/data/P_Mfchif_r.dat.gz', P_Mfchif_r)
np.savetxt(out_dir+'/data/P_Mfchif_imr.dat.gz', P_Mfchif_imr)
np.savetxt(out_dir+'/data/dMfbyMf_vec.dat.gz', dMfbyMf_vec)
np.savetxt(out_dir+'/data/dchifbychif_vec.dat.gz', dchifbychif_vec)
np.savetxt(out_dir+'/data/P_dMfbyMf_dchifbychif.dat.gz', P_dMfbyMf_dchifbychif)
np.savetxt(out_dir+'/data/P_dMfbyMf.dat.gz', P_dMfbyMf)
np.savetxt(out_dir+'/data/P_dchifbychif.dat.gz', P_dchifbychif)
np.savetxt(out_dir+'/data/GR_confidence.txt', [gr_conf_level])

#--------------------------------------------------------------------------------------
# Mean values and 90% intervals for final mass and spin from the combined IMR posterior
#--------------------------------------------------------------------------------------
P_Mf = np.sum(P_Mfchif_imr, axis=0) * dchif_intp
P_chif = np.sum(P_Mfchif_imr, axis=1) * dMf_intp

conf_Mf = confidence(P_Mf)
s_Mf = conf_Mf.height_from_level(0.9)

conf_chif = confidence(P_chif)
s_chif = conf_chif.height_from_level(0.9)

left_Mf = min(Mf_intp[np.where(P_Mf>=s_Mf)[0]])
right_Mf = max(Mf_intp[np.where(P_Mf>=s_Mf)[0]])

left_chif = min(chif_intp[np.where(P_chif>=s_chif)[0]])
right_chif = max(chif_intp[np.where(P_chif>=s_chif)[0]])

mean_Mf = np.average(Mf_intp, weights=P_Mf)
mean_chif = np.average(chif_intp, weights=P_chif)

np.savetxt(out_dir+'/data/Mf_values.txt',[mean_Mf, left_Mf-mean_Mf, right_Mf-mean_Mf])
np.savetxt(out_dir+'/data/chif_values.txt',[mean_chif, left_chif-mean_chif, right_chif-mean_chif])
