#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

seobnr_raw_pos = np.genfromtxt('/home/abhirup/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DS_nonprecspin_Healy2014_2016-01-01_fhigh_insp132Hz_flow_ring_132Hz_final_run/lalinf_imr/posterior_samples.dat', dtype=None, names=True)

imrpp_raw_pos = np.genfromtxt('/home/abhirup/public_html/imrtestgr/ER8/G184098/IMRPhenomPv2_nonprecspin_Healy2014_2016-01-01_fhigh_insp132Hz_flow_ring_132Hz_final_run/lalinf_imr/posterior_samples.dat', dtype=None, names=True)

seobnr_pos = np.loadtxt('/home/abhirup/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DS_nonprecspin_Healy2014_2016-01-01_fhigh_insp132Hz_flow_ring_132Hz_final_run/data/P_Mfchif_imr.dat.gz')

(seobnr_Mf_edg, seobnr_chif_edg) = np.loadtxt('/home/abhirup/public_html/imrtestgr/ER8/G184098/SEOBNRv2_ROM_DS_nonprecspin_Healy2014_2016-01-01_fhigh_insp132Hz_flow_ring_132Hz_final_run/data/Mfchif.dat.gz')

imrpp_pos = np.loadtxt('/home/abhirup/public_html/imrtestgr/ER8/G184098/IMRPhenomPv2_nonprecspin_Healy2014_2016-01-01_fhigh_insp132Hz_flow_ring_132Hz_final_run/data/P_Mfchif_imr.dat.gz')

(imrpp_Mf_edg, imrpp_chif_edg) = np.loadtxt('/home/abhirup/public_html/imrtestgr/ER8/G184098/IMRPhenomPv2_nonprecspin_Healy2014_2016-01-01_fhigh_insp132Hz_flow_ring_132Hz_final_run/data/Mfchif.dat.gz')

seobnr_Mf_cen = (seobnr_Mf_edg[1:] + seobnr_Mf_edg[:-1])/2.
seobnr_chif_cen = (seobnr_chif_edg[1:] + seobnr_chif_edg[:-1])/2.
seobnr_dMf = np.mean(np.diff(seobnr_Mf_cen))
seobnr_Mf_pos = np.sum(seobnr_pos, axis=0)
seobnr_Mf_integral = np.sum(seobnr_Mf_pos) * seobnr_dMf
seobnr_Mf_norm_pos = seobnr_Mf_pos / seobnr_Mf_integral
seobnr_dchif = np.mean(np.diff(seobnr_chif_cen))
seobnr_chif_pos = np.sum(seobnr_pos, axis=1)
seobnr_chif_integral = np.sum(seobnr_chif_pos) * seobnr_dchif
seobnr_chif_norm_pos = seobnr_chif_pos / seobnr_chif_integral


seobnr_Mf_mean = np.average(seobnr_Mf_cen, weights=np.sum(seobnr_pos, axis=0))
seobnr_chif_mean = np.average(seobnr_chif_cen, weights=np.sum(seobnr_pos, axis=1))

seobnr_Mf_cumarr = np.cumsum(seobnr_Mf_pos)/np.sum(seobnr_Mf_pos)
seobnr_inverse_Mf_cuminterp = interpolate.interp1d(seobnr_Mf_cumarr, seobnr_Mf_cen)
seobnr_Mf_median = seobnr_inverse_Mf_cuminterp(0.5)
seobnr_chif_cumarr = np.cumsum(seobnr_chif_pos)/np.sum(seobnr_chif_pos)
seobnr_inverse_chif_cuminterp = interpolate.interp1d(seobnr_chif_cumarr, seobnr_chif_cen)
seobnr_chif_median = seobnr_inverse_chif_cuminterp(0.5)


plt.figure()
(counts, bins, patches) = plt.hist(seobnr_raw_pos['mf'], bins=100, normed=True, histtype='step', color='r', label='Initial')
plt.plot(seobnr_Mf_cen, seobnr_Mf_norm_pos, 'g-', label='Final')
plt.xlim(55, 80)
plt.legend(loc='best')
plt.savefig('../plots/seobnr_Mf_prior_comp.png')

plt.figure()
(counts, bins, patches) = plt.hist(seobnr_raw_pos['af'], bins=100, normed=True, histtype='step', color='r', label='Initial')
plt.plot(seobnr_chif_cen, seobnr_chif_norm_pos, 'g-', label='Final')
plt.xlim(0.45, 0.85)
plt.legend(loc='best')
plt.savefig('../plots/seobnr_chif_prior_comp.png')

imrpp_Mf_cen = (imrpp_Mf_edg[1:] + imrpp_Mf_edg[:-1])/2.
imrpp_chif_cen = (imrpp_chif_edg[1:] + imrpp_chif_edg[:-1])/2.
imrpp_dMf = np.mean(np.diff(imrpp_Mf_cen))
imrpp_Mf_pos = np.sum(imrpp_pos, axis=0)
imrpp_Mf_integral = np.sum(imrpp_Mf_pos) * imrpp_dMf
imrpp_Mf_norm_pos = imrpp_Mf_pos / imrpp_Mf_integral
imrpp_dchif = np.mean(np.diff(imrpp_chif_cen))
imrpp_chif_pos = np.sum(imrpp_pos, axis=1)
imrpp_chif_integral = np.sum(imrpp_chif_pos) * imrpp_dchif
imrpp_chif_norm_pos = imrpp_chif_pos / imrpp_chif_integral


imrpp_Mf_mean = np.average(imrpp_Mf_cen, weights=np.sum(imrpp_pos, axis=0))
imrpp_chif_mean = np.average(imrpp_chif_cen, weights=np.sum(imrpp_pos, axis=1))

imrpp_Mf_cumarr = np.cumsum(imrpp_Mf_pos)/np.sum(imrpp_Mf_pos)
imrpp_inverse_Mf_cuminterp = interpolate.interp1d(imrpp_Mf_cumarr, imrpp_Mf_cen)
imrpp_Mf_median = imrpp_inverse_Mf_cuminterp(0.5)
imrpp_chif_cumarr = np.cumsum(imrpp_chif_pos)/np.sum(imrpp_chif_pos)
imrpp_inverse_chif_cuminterp = interpolate.interp1d(imrpp_chif_cumarr, imrpp_chif_cen)
imrpp_chif_median = imrpp_inverse_chif_cuminterp(0.5)



plt.figure()
(counts, bins, patches) = plt.hist(imrpp_raw_pos['mf'], bins=100, normed=True, histtype='step', color='r', label='Initial')
plt.plot(imrpp_Mf_cen, imrpp_Mf_norm_pos, 'g-', label='Final')
plt.xlim(55, 80)
plt.legend(loc='best')
plt.savefig('../plots/imrpp_Mf_prior_comp.png')

plt.figure()
(counts, bins, patches) = plt.hist(imrpp_raw_pos['af'], bins=100, normed=True, histtype='step', color='r', label='Initial')
plt.plot(imrpp_chif_cen, imrpp_chif_norm_pos, 'g-', label='Final')
plt.xlim(0.45, 0.85)
plt.legend(loc='best')
plt.savefig('../plots/imrpp_chif_prior_comp.png')

