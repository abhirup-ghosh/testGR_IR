import matplotlib as mpl
mpl.use('Agg')
import os, sys, numpy as np
import matplotlib.pyplot as plt
import bayesian as ba
import standard_gwtransf as gw
import imrtgrutils_final as tgr
import scipy.ndimage.filters as filter

def gf(P):
        return filter.gaussian_filter(P, sigma=4.0)

indir_root = '/home/abhirup/Documents/Work/testGR_IR/runs/td_likelihood/20170925_semi-final_runs/20170925_GW150914'
post_loc_1 = indir_root + '/SEOBNRv4_opt_m1_50_m2_50_a1z_0.00_a2z_0.00_dL_500_flow_30_srate_2048_CE_t0_ring_3ms_QNM'
post_loc_2 = indir_root + '/SEOBNRv4_opt_m1_50_m2_50_a1z_0.00_a2z_0.00_dL_500_flow_30_srate_2048_CE_t0_ring_5ms_QNM'
post_loc_3 = indir_root + '/SEOBNRv4_opt_m1_50_m2_50_a1z_0.00_a2z_0.00_dL_500_flow_30_srate_2048_CE_t0_ring_10ms_QNM'

post_loc_list = [post_loc_1, post_loc_2, post_loc_3]
method = 'samples'
QNM = True
color = ['r', 'k', 'g']

plt.figure(figsize=(16,5))
for (i,post_loc) in enumerate(post_loc_list):
  if method == 'imrtgr':
	Mf_bins,af_bins = np.loadtxt(post_loc+'/Mfchif.dat.gz')
	P_Mfaf_i = np.loadtxt(post_loc+'/P_Mfchif_i.dat.gz')
	P_Mfaf_r = np.loadtxt(post_loc+'/P_Mfchif_r.dat.gz')

	eps_vec = np.loadtxt(post_loc_FD+'/dMfbyMf_vec.dat.gz')
	sig_vec = np.loadtxt(post_loc_FD+'/dchifbychif_vec.dat.gz')
	P_epssig = np.loadtxt(post_loc+'/P_dMfbyMf_dchifbychif.dat.gz')

  if method == 'samples':
	insp_data = np.genfromtxt(post_loc+'/pe/data/inspiral_samples_postburnin.dat', names=True, dtype=None)
	ring_data = np.genfromtxt(post_loc+'/pe/data/ringdown_samples_postburnin.dat', names=True, dtype=None)
	Mf_i, af_i = insp_data['mf'], insp_data['af']
	Mf_r, af_r = ring_data['mf'], ring_data['af']
	Mf_bins = np.linspace(min(min(Mf_i), min(Mf_r)), max(max(Mf_i), max(Mf_r)), 201)
	af_bins = np.linspace(min(min(af_i), min(af_r)), max(max(af_i), max(af_r)), 201)
	P_Mfaf_i, Mf_bins, af_bins = np.histogram2d(Mf_i, af_i, bins=(Mf_bins, af_bins), normed=True)
	P_Mfaf_r, Mf_bins, af_bins = np.histogram2d(Mf_r, af_r, bins=(Mf_bins, af_bins), normed=True)
	P_Mfaf_i = P_Mfaf_i.T
	P_Mfaf_r = P_Mfaf_r.T
	if QNM == True:
        	tau, f_qnm = ring_data['tau'], ring_data['f_qnm']
        	tau_bins = np.linspace(min(tau), max(tau), 201)
        	f_qnm_bins = np.linspace(min(f_qnm), max(f_qnm), 201)
       		P_tauf_qnm, tau_bins, f_qnm_bins = np.histogram2d(tau, f_qnm, bins=(tau_bins, f_qnm_bins), normed=True)
        	P_tauf_qnm = P_tauf_qnm.T
		s1_tf, s2_tf = ba.nsigma_value(P_tauf_qnm, 0.5), ba.nsigma_value(P_tauf_qnm, 0.9)

	eps_vec = np.linspace(-1, 1, 201)
	sig_vec = np.linspace(-1, 1, 201)
	P_epssig, eps_vec, sig_vec = np.histogram2d((Mf_i-Mf_r)/(Mf_i+Mf_r), (af_i-af_r)/(af_i+af_r), bins=(eps_vec, sig_vec), normed=True)
	P_epssig = P_epssig.T

  s1_i, s2_i = ba.nsigma_value(P_Mfaf_i, 0.5), ba.nsigma_value(P_Mfaf_i, 0.9)
  s1_r, s2_r = ba.nsigma_value(P_Mfaf_r, 0.5), ba.nsigma_value(P_Mfaf_r, 0.9)
  s1, s2 = ba.nsigma_value(P_epssig, 0.5), ba.nsigma_value(P_epssig, 0.9)

  plt.subplot(132)
  plt.contour(Mf_bins[:-1],af_bins[:-1],P_Mfaf_i, levels=(s1_i, s2_i), colors=color[i])
  plt.contour(Mf_bins[:-1],af_bins[:-1],P_Mfaf_r, levels=(s1_r, s2_r), ls=('--','--'), colors=color[i])
  plt.axvline(x=95.15)
  plt.axhline(y=0.69)
  plt.xlabel('$M_f$')
  plt.ylabel('$a_f$')
  plt.xlim([90, 130])
  plt.ylim([0.5,0.9])
  if QNM == True:
        plt.subplot(131)
        plt.contour(tau_bins[:-1], f_qnm_bins[:-1], P_tauf_qnm, levels=(s1_tf, s2_tf), colors=color[i])
  	plt.xlim([0.004, 0.01])
  	plt.ylim([150,175])
	plt.xlabel('$\\tau$')
	plt.ylabel('$f_{QNM}$')
  plt.subplot(133)
  plt.contour(eps_vec[:-1], sig_vec[:-1], P_epssig, levels=(s1, s2), colors=color[i])
  plt.xlabel('$\Delta M_f/M_f$')
  plt.ylabel('$\Delta a_f/a_f$')
  plt.xlim([-0.2, 0.1])
  plt.ylim([-0.15, 0.1])
  plt.hold(True)
  
  plt.tight_layout()


plt.savefig('../plots/td_likelihood/TD_QNM_different_t0_CE.png', dpi=300)

