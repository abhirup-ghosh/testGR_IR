import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import commands, os
import utility_codes as uc
import imrtestgr as tgr
import nr_fits as nr
import pickle, gzip
import glob
import extract_from_lalinferencenest_hdf5 as exh5

inj_list = np.loadtxt('common_injection_list.txt')
fit_formula = 'nonprecspin_Healy2014'
N_bins = 201

method = 2

prior_Mfaf_file = '/home/abhirup/opt/lalsuite_lalinference_o2/share/lalinference/imrtgr_prior_data/Prior_Mfaf_nonprec_Healy2014_M1-500_aligned.pklz'



if method == 1:
  B_GR = []
  B_modGR = []
  B_modGR_var = []
  for inj_type in ['GR', 'modGR_a2_20', 'modGR_variable_a2']:   
    for inj in inj_list[0:100]:

      # selecting the appropriate posterior samples location	
      if inj_type == 'GR':
      	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume/imrtestgr_nonprecspin_Healy2014_bins_401_bins_-2.0_2.0_w_priorcorr/injection_%d'%inj
      elif inj_type == 'modGR_a2_20':
	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/modGR_a2_20/imrtestgr_nonprecspin_Healy2014_bins_401_-2.0_2.0_w_priorcorr/IHES_%04d'%inj
      elif inj_type == 'modGR_variable_a2':
	post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/modGR_variable_a2/imrtestgr_nonprecspin_Healy2014_bins_401_-2.0_2.0_w_priorcorr/IHES_%04d'%inj

      
      # read the text files containing the optimal snr
      insp_summary_stats_loc = '%s/lalinf_insp/summary_statistics.dat'%post_loc
      ring_summary_stats_loc = '%s/lalinf_ring/summary_statistics.dat'%post_loc

      insp_snr, insp_h1_snr, insp_l1_snr, insp_v1_snr, insp_m1_rec, insp_m2_rec, insp_chi1_rec, insp_chi2_rec, insp_Mf_rec, insp_chif_rec = uc.optimal_snr_module(insp_summary_stats_loc)
      ring_snr, ring_h1_snr, ring_l1_snr, ring_v1_snr, ring_m1_rec, ring_m2_rec, ring_chi1_rec, ring_chi2_rec, ring_Mf_rec, ring_chif_rec = uc.optimal_snr_module(ring_summary_stats_loc)


      if insp_snr > 8 and ring_snr > 8:

	# Reading data from each posterior location
	data_i = np.genfromtxt('%s/lalinf_insp/posterior_samples.dat'%post_loc, names=True, dtype=None)
	data_r = np.genfromtxt('%s/lalinf_ring/posterior_samples.dat'%post_loc, names=True, dtype=None)

	Mf_i, af_i = data_i['mf'], data_i['af']
	Mf_r, af_r = data_r['mf'], data_r['af']

	# fixing the bounds of the histogram for binning the mf, af samples
	Mf_min = min(min(Mf_i), min(Mf_r))
	Mf_max = max(max(Mf_i), max(Mf_r))
	af_min = min(min(af_i), min(af_r))
        af_max = max(max(af_i), max(af_r))
	
	Mf_bins = np.linspace(Mf_min, Mf_max, N_bins)
	af_bins = np.linspace(af_min, af_max, N_bins)

	dMf = np.mean(np.diff(Mf_bins))
	daf = np.mean(np.diff(af_bins))

	Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
	af_intp = (af_bins[:-1] + af_bins[1:])/2.

	# binning mf, af samples over common histogram
	P_Mfaf_i, Mf_bins, af_bins = np.histogram2d(Mf_i, af_i, bins=(Mf_bins, af_bins), normed=True)
	P_Mfaf_r, Mf_bins, af_bins = np.histogram2d(Mf_r, af_r, bins=(Mf_bins, af_bins), normed=True)

	f = gzip.open(prior_Mfaf_file,'rb')
	P_Mfaf_pr_interp_obj = pickle.load(f)
	P_Mfaf_pr = P_Mfaf_pr_interp_obj(Mf_intp, af_intp)

	P_Mfaf_i = P_Mfaf_i.T
	P_Mfaf_r = P_Mfaf_r.T

	# normalising the posteriors
	P_Mfaf_i /= np.sum(P_Mfaf_i) * dMf * daf
	P_Mfaf_r /= np.sum(P_Mfaf_r) * dMf * daf
	P_Mfaf_pr /= np.sum(P_Mfaf_pr) * dMf * daf

	B_integrand = P_Mfaf_i * P_Mfaf_r / P_Mfaf_pr
	B_integrand[np.isnan(B_integrand)] = 0.
	B_integrand[np.isinf(B_integrand)] = 0.

	B_post_loc =  np.sum(B_integrand) * dMf * daf 

	print '... Bayes factor for %s injection %d: %f'%(inj_type, inj, B_post_loc)

	# sorting the Bayes factor in the appropriate array for GR or modGR values
	if inj_type == 'GR':
	  B_GR = np.append(B_GR, B_post_loc)
	elif inj_type == 'modGR_a2_20':
          B_modGR = np.append(B_modGR, B_post_loc)
	elif inj_type == 'modGR_variable_a2':
	  B_modGR_var = np.append(B_modGR_var, B_post_loc)

  print B_GR
  print B_modGR
  print B_modGR_var

  plt.figure()
  for idx in range(len(B_GR)):
	plt.plot(idx + 1, np.log(B_GR[idx]), marker='o', markerfacecolor='none', markeredgecolor='r')
	plt.plot(idx + 1, np.log(np.prod(B_GR[0:idx])), marker='x', color='r')
	plt.plot(idx + 1, np.log(B_modGR[idx]), marker='o', markerfacecolor='none', color='k')
        plt.plot(idx + 1, np.log(np.prod(B_modGR[0:idx])), marker='x', color='k')
	plt.plot(idx + 1, np.log(B_modGR_var[idx]), marker='o', markerfacecolor='none', color='g')
        plt.plot(idx + 1, np.log(np.prod(B_modGR_var[0:idx])), marker='x', color='g')
	plt.hold(True)
  plt.legend(['GR combined', 'GR individual', 'modGR combined', 'modGR individual', 'var modGR combined', 'var modGR individual'], loc='best')
  plt.ylabel('log B')
  plt.xlabel('no of injections')
  plt.savefig('../plots/bayes_factor_priornormalised.png')


# alternate method
if method == 2:
  O_GR = []
  O_modGR = []
  O_modGR_var = []
  for inj_type in ['GR', 'modGR_a2_20', 'modGR_variable_a2']:
    for inj in inj_list[0:100]:

      # selecting the appropriate posterior samples location
      if inj_type == 'GR':
        post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-01-16/uniform_compmass_spins_comoving_volume'

	# read the text files containing the optimal snr
        insp_summary_stats_loc = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_bins_-2.0_2.0_w_priorcorr/injection_%d/lalinf_insp/summary_statistics.dat'%(post_loc, inj)
        ring_summary_stats_loc = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_bins_-2.0_2.0_w_priorcorr/injection_%d/lalinf_ring/summary_statistics.dat'%(post_loc, inj)

        insp_snr, insp_h1_snr, insp_l1_snr, insp_v1_snr, insp_m1_rec, insp_m2_rec, insp_chi1_rec, insp_chi2_rec, insp_Mf_rec, insp_chif_rec = uc.optimal_snr_module(insp_summary_stats_loc)
        ring_snr, ring_h1_snr, ring_l1_snr, ring_v1_snr, ring_m1_rec, ring_m2_rec, ring_chi1_rec, ring_chi2_rec, ring_Mf_rec, ring_chif_rec = uc.optimal_snr_module(ring_summary_stats_loc)

	if insp_snr > 8 and ring_snr > 8:

	  B_imr = np.genfromtxt(glob.glob('%s/injection_%d/IMR/engine/*.dat_B.txt'%(post_loc, inj))[0])[0:1]
	  B_insp = np.genfromtxt(glob.glob('%s/injection_%d/inspiral/engine/*.dat_B.txt'%(post_loc, inj))[0])[0:1]
	  B_ring = np.genfromtxt(glob.glob('%s/injection_%d/post-inspiral/engine/*.dat_B.txt'%(post_loc, inj))[0])[0:1]
	  O_GR_inj = np.exp(B_insp - B_imr) + np.exp(B_ring - B_imr)
	  O_GR = np.append(O_GR, O_GR_inj)

      elif inj_type == 'modGR_a2_20':
        post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/modGR_a2_20'

	# read the text files containing the optimal snr
        insp_summary_stats_loc = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_-2.0_2.0_w_priorcorr/IHES_%04d/lalinf_insp/summary_statistics.dat'%(post_loc, inj)
        ring_summary_stats_loc = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_-2.0_2.0_w_priorcorr/IHES_%04d/lalinf_ring/summary_statistics.dat'%(post_loc, inj)

        insp_snr, insp_h1_snr, insp_l1_snr, insp_v1_snr, insp_m1_rec, insp_m2_rec, insp_chi1_rec, insp_chi2_rec, insp_Mf_rec, insp_chif_rec = uc.optimal_snr_module(insp_summary_stats_loc)
        ring_snr, ring_h1_snr, ring_l1_snr, ring_v1_snr, ring_m1_rec, ring_m2_rec, ring_chi1_rec, ring_chi2_rec, ring_Mf_rec, ring_chif_rec = uc.optimal_snr_module(ring_summary_stats_loc)

	if insp_snr > 8 and ring_snr > 8:

	  B_imr = exh5.get_metadata(glob.glob('%s/IHES_%04d/IMR/*/*/engine/*.hdf5'%(post_loc, inj))[0])['log_bayes_factor']
	  B_insp = exh5.get_metadata(glob.glob('%s/IHES_%04d/inspiral/*/*/engine/*.hdf5'%(post_loc, inj))[0])['log_bayes_factor']
	  B_ring = exh5.get_metadata(glob.glob('%s/IHES_%04d/post-inspiral/*/*/engine/*.hdf5'%(post_loc, inj))[0])['log_bayes_factor']
          O_modGR_inj = np.exp(B_insp - B_imr) + np.exp(B_ring - B_imr)
	  O_modGR = np.append(O_modGR, O_modGR_inj)
	  

      elif inj_type == 'modGR_variable_a2':
        post_loc = '/home/abhirup/Documents/Work/testGR_IR/runs/simulations_modGR/gaussian_noise/LALAdLIGO/SEOBNRv2_ROM_DoubleSpinthreePointFivePN/2016-11-30_3det/modGR_variable_a2'

	# read the text files containing the optimal snr
        insp_summary_stats_loc = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_-2.0_2.0_w_priorcorr/IHES_%04d/lalinf_insp/summary_statistics.dat'%(post_loc, inj)
        ring_summary_stats_loc = '%s/imrtestgr_nonprecspin_Healy2014_bins_401_-2.0_2.0_w_priorcorr/IHES_%04d/lalinf_ring/summary_statistics.dat'%(post_loc, inj)

        insp_snr, insp_h1_snr, insp_l1_snr, insp_v1_snr, insp_m1_rec, insp_m2_rec, insp_chi1_rec, insp_chi2_rec, insp_Mf_rec, insp_chif_rec = uc.optimal_snr_module(insp_summary_stats_loc)
        ring_snr, ring_h1_snr, ring_l1_snr, ring_v1_snr, ring_m1_rec, ring_m2_rec, ring_chi1_rec, ring_chi2_rec, ring_Mf_rec, ring_chif_rec = uc.optimal_snr_module(ring_summary_stats_loc)

        if insp_snr > 8 and ring_snr > 8:

          B_imr = exh5.get_metadata(glob.glob('%s/IHES_%04d/IMR/*/*/engine/*.hdf5'%(post_loc, inj))[0])['log_bayes_factor']
          B_insp = exh5.get_metadata(glob.glob('%s/IHES_%04d/inspiral/*/*/engine/*.hdf5'%(post_loc, inj))[0])['log_bayes_factor']
          B_ring = exh5.get_metadata(glob.glob('%s/IHES_%04d/post-inspiral/*/*/engine/*.hdf5'%(post_loc, inj))[0])['log_bayes_factor']
	  O_modGR_var_inj = np.exp(B_insp - B_imr) + np.exp(B_ring - B_imr)
          O_modGR_var = np.append(O_modGR_var, O_modGR_var_inj)

  idx1, = np.where(O_GR == 0.)
  idx2, = np.where(O_modGR == 0.)
  idx3, = np.where(O_modGR_var == 0.)
  idx = np.append(np.append(idx1, idx2), idx3)

  O_GR = np.delete(O_GR, idx)
  O_modGR = np.delete(O_modGR, idx)
  O_modGR_var = np.delete(O_modGR_var, idx)

  O_GR = np.log(O_GR)
  O_modGR = np.log(O_modGR)
  O_modGR_var = np.log(O_modGR_var)

  plt.figure()
  plt.plot(np.arange(len(O_GR)) + 1, np.cumsum(O_GR), marker='x', color='r', ls='None')
  plt.plot(np.arange(len(O_modGR)) + 1, np.cumsum(O_modGR), marker='x', color='k', ls='None')
  plt.plot(np.arange(len(O_modGR_var)) + 1, np.cumsum(O_modGR_var), marker='x', color='c', ls='None')
  plt.legend(['GR combined', 'GR individual', 'modGR combined', 'modGR individual', 'var modGR combined', 'var modGR individual'], loc='best')
  plt.ylabel('log B')
  plt.xlabel('no of injections')
  plt.savefig('../plots/bayes_factor_method2_combined.png')

  plt.figure()
  plt.plot(np.arange(len(O_GR)) + 1, O_GR, marker='x', color='r', ls='None')
  plt.plot(np.arange(len(O_modGR)) + 1, O_modGR, marker='x', color='k', ls='None')
  plt.plot(np.arange(len(O_modGR_var)) + 1, O_modGR_var, marker='x', color='c', ls='None')
  plt.legend(['GR', 'modGR', 'var modGR'], loc='best')
  plt.ylabel('log B')
  plt.xlabel('no of injections')
  plt.savefig('../plots/bayes_factor_method2_casewise.png')
