import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import bayesian as ba
import scipy.ndimage.filters as filter

""" gaussian filter of histogram """
def gf(P):
        return filter.gaussian_filter(P, sigma=4.0)

data_i = np.genfromtxt('inspiral_samples_postburnin.dat', names=True, dtype=None)
mf_i, af_i, lnprob_i = data_i['mf'], data_i['af'], data_i['lnprob']
data_r = np.genfromtxt('ringdown_samples_postburnin.dat', names=True, dtype=None)
mf_r, af_r, lnprob_r = data_r['mf'], data_r['af'], data_r['lnprob']

mf_bins = np.linspace(min(min(mf_i), min(mf_r)), max(max(mf_i), max(mf_r)), 401)
af_bins = np.linspace(min(min(af_i), min(af_r)), max(max(af_i), max(af_r)), 401)

P_mfaf_i, mf_bins, af_bins = np.histogram2d(mf_i, af_i, bins=(mf_bins, af_bins))
P_mfaf_i = P_mfaf_i.T
P_mfaf_r, mf_bins, af_bins = np.histogram2d(mf_r, af_r, bins=(mf_bins, af_bins))
P_mfaf_r = P_mfaf_r.T

s1_i, s2_i = ba.nsigma_value(P_mfaf_i, 0.68), ba.nsigma_value(P_mfaf_i, 0.95)
s1_r, s2_r = ba.nsigma_value(P_mfaf_r, 0.68), ba.nsigma_value(P_mfaf_r, 0.95)

print len(mf_bins), len(af_bins), np.shape(P_mfaf_i)

plt.figure(figsize=(10,6))
plt.subplot(131)
plt.scatter(mf_i, af_i, c=lnprob_i, lw=0)
plt.contour(mf_bins[:-1], af_bins[:-1], gf(P_mfaf_i), levels=(s1_i, s2_i), colors='orange')
plt.subplot(132)
plt.scatter(mf_r, af_r, c=lnprob_r, lw=0)
plt.contour(mf_bins[:-1], af_bins[:-1], gf(P_mfaf_r), levels=(s1_r, s2_r), colors='r')
plt.subplot(133)
plt.contour(mf_bins[:-1], af_bins[:-1], gf(P_mfaf_i), levels=(s1_i, s2_i), colors='orange')
plt.contour(mf_bins[:-1], af_bins[:-1], gf(P_mfaf_r), levels=(s1_r, s2_r), colors='r')
plt.savefig('m1m2.png')

