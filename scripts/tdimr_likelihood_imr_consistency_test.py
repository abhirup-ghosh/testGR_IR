#import matplotlib as mpl
#mpl.use('Agg')
import os, sys, numpy as np
import matplotlib.pyplot as plt
import bayesian as ba

post_loc_TD = '/Users/abhirupghosh/Documents/Work/testGR_IR/runs/td_likelihood/20170925_semi-final_runs/20170925_GW150914/SEOBNRv4_opt_m1_36_m2_29_a1z_0.00_a2z_0.00_dL_500_flow_30_srate_2048_AdLIGO_ZDHP_t0_ring_0ms_lalsimulation/pe/data'
post_loc_FD = '/Users/abhirupghosh/Documents/Work/testGR_IR/runs/td_likelihood/20170925_semi-final_runs/20170925_GW150914/SEOBNRv4_opt_m1_36_m2_29_a1z_0.00_a2z_0.00_dL_500_flow_30_srate_2048_ZERO_DET_high_P_t0_ring_0ms_FD'

insp_data_TD = np.genfromtxt(post_loc_TD+'/inspiral_samples_postburnin.dat',dtype=None,names=True)
ring_data_TD = np.genfromtxt(post_loc_TD+'/ringdown_samples_postburnin.dat',dtype=None,names=True)
insp_data_FD = np.genfromtxt(post_loc_FD+'/inspiral_H1_132Hz_pinparams/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat',dtype=None,names=True)
ring_data_FD = np.genfromtxt(post_loc_FD+'/post-inspiral_H1_132Hz_pinparams/lalinferencenest/SEOBNRv4_ROMpseudoFourPN/1126285216.000000000-0/H1/posterior_samples.dat',dtype=None,names=True)

mf_i_TD, af_i_TD = insp_data_TD['mf'], insp_data_TD['af']
mf_r_TD, af_r_TD = ring_data_TD['mf'], ring_data_TD['af']
mf_i_FD, af_i_FD = insp_data_FD['mf'], insp_data_FD['af']
mf_r_FD, af_r_FD = ring_data_FD['mf'], ring_data_FD['af']

print len(mf_i_FD), len(af_i_FD)
print len(mf_r_FD), len(af_r_FD)

mf_bins = np.linspace(40, 80, 101)
af_bins = np.linspace(0, 1, 101)

P_mfaf_i_TD, mf_bins, af_bins = np.histogram2d(mf_i_TD, af_i_TD, bins=(mf_bins, af_bins), normed=True)
P_mfaf_r_TD, mf_bins, af_bins = np.histogram2d(mf_r_TD, af_r_TD, bins=(mf_bins, af_bins), normed=True)
P_mfaf_i_FD, mf_bins, af_bins = np.histogram2d(mf_i_FD, af_i_FD, bins=(mf_bins, af_bins), normed=True)
P_mfaf_r_FD, mf_bins, af_bins = np.histogram2d(mf_r_FD, af_r_FD, bins=(mf_bins, af_bins), normed=True)

P_mfaf_i_TD = P_mfaf_i_TD.T
P_mfaf_r_TD = P_mfaf_r_TD.T
P_mfaf_i_FD = P_mfaf_i_FD.T
P_mfaf_r_FD = P_mfaf_r_FD.T

s1_i_TD, s2_i_TD = ba.nsigma_value(P_mfaf_i_TD, 0.5), ba.nsigma_value(P_mfaf_i_TD, 0.9)
s1_r_TD, s2_r_TD = ba.nsigma_value(P_mfaf_r_TD, 0.5), ba.nsigma_value(P_mfaf_r_TD, 0.9)
s1_i_FD, s2_i_FD = ba.nsigma_value(P_mfaf_i_FD, 0.5), ba.nsigma_value(P_mfaf_i_FD, 0.9)
s1_r_FD, s2_r_FD = ba.nsigma_value(P_mfaf_r_FD, 0.5), ba.nsigma_value(P_mfaf_r_FD, 0.9)


eps_bins = np.linspace(-1, 1, 51)
sig_bins = np.linspace(-1, 1, 51)

P_epssig_TD, eps_bins, sig_bins = np.histogram2d(mf_i_TD - mf_r_TD, af_i_TD - af_i_TD, bins=(eps_bins, sig_bins), normed=True)
P_epssig_FD, eps_bins, sig_bins = np.histogram2d(mf_i_FD[-10000:] - mf_r_FD[-10000:], af_i_FD[-10000:] - af_i_FD[-10000:], bins=(eps_bins, sig_bins), normed=True)

P_epssig_TD = P_epssig_TD.T
P_epssig_FD = P_epssig_FD.T

s1_TD, s2_TD = ba.nsigma_value(P_epssig_TD, 0.5), ba.nsigma_value(P_epssig_TD, 0.9)
s1_FD, s2_FD = ba.nsigma_value(P_epssig_FD, 0.5), ba.nsigma_value(P_epssig_FD, 0.9)


plt.figure(figsize=(10,5))
plt.subplot(121)
plt.contour(mf_bins[:-1], af_bins[:-1], P_mfaf_i_TD, levels=(s1_i_TD, s2_i_TD), colors='r')
plt.contour(mf_bins[:-1], af_bins[:-1], P_mfaf_r_TD, levels=(s1_r_TD, s2_r_TD), colors='orange')
plt.contour(mf_bins[:-1], af_bins[:-1], P_mfaf_i_FD, levels=(s1_i_FD, s2_i_FD),colors='b')
plt.contour(mf_bins[:-1], af_bins[:-1], P_mfaf_r_FD, levels=(s1_r_FD, s2_r_FD),colors='g')
plt.xlim([55, 70])
plt.ylim([0.6,0.75])
plt.subplot(122)
plt.contour(eps_bins[:-1], sig_bins[:-1], P_epssig_TD, levels=(s1_TD, s2_TD), colors='r')
plt.contour(eps_bins[:-1], sig_bins[:-1], P_epssig_FD, levels=(s1_FD, s2_FD), colors='k')
plt.xlim([-1,1])
plt.ylim([-0.01,0.01])
plt.show()
