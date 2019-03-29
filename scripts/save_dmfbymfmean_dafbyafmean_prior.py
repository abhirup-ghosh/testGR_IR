""" 
Compute the prior distribution in the fractional mass and spin parameters (dMfbyMf, dafbyaf) 
corresponding to a uniform prior in final and spin of the final black hole. 

Abhirup Ghosh, 2016-06-24 

$Id:$
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import tgrplotsettings, time
import os
sys.path.insert(1, os.path.join(os.path.expanduser('~'), 'src/lalsuite/lalinference/python/lalinference/imrtgr'))
import imrtgrutils as tgr
import scipy
import scipy.signal as ss
from scipy import interpolate
import pickle, gzip
from optparse import OptionParser

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

if __name__ == '__main__':

  start_time = time.time()

  # read inputs from command line 
  parser = OptionParser()
  parser.add_option("--fit-formula", dest="fit_formula", help="fitting formula for the mass/spin of the final BH [options: 'nospin_Pan2011', 'nonprecspin_Healy2014', 'bbh_average_fits_precessing']", default="nonprecspin_Healy2014")
  parser.add_option("--spin_angle_dist", dest="spin_angle_dist", help="[options: 'aligned', 'isotropic']", default='aligned')
  parser.add_option("--uniform-in", dest="uniform_in", 
		  help="tag to define the random samples uniform in (m1,m2,a1,a2) or (Mf, af) [options=:'m1m2a1a2', 'Mfaf']", default='m1m2a1a2')
  parser.add_option("--comp-mass-min", dest="comp_mass_min",
                  help="minimum value of the component mass prior [M_sun]", default=1)
  parser.add_option("--comp-mass-max", dest="comp_mass_max",
                  help="maximum value of the component mass prior [M_sun]", default=300)
  parser.add_option("--comp-spin-min", dest="comp_spin_min",
                  help="minimum value of the dimensionless spin prior", default=0)
  parser.add_option("--comp-spin-max", dest="comp_spin_max",
                  help="maximum value of the dimensionless spin prior", default=1) 
  parser.add_option("--Mf-min", dest="Mf_min",
                  help="minimum value of the Mf prior [M_sun]", default=1)
  parser.add_option("--Mf-max", dest="Mf_max",
                  help="maximum value of the Mf prior [M_sun]", default=1000)
  parser.add_option("--af-min", dest="af_min",
                  help="minimum value of the af prior", default=0)
  parser.add_option("--af-max", dest="af_max",
                  help="maximum value of the af prior", default=1)
  parser.add_option("--dMfbyMf-lim", dest="dMfbyMf_lim",
                  help="the prior distribution will be calculated from -dMfbyMf_lim to dMfbyMf_lim", default=1)
  parser.add_option("--dafbyaf-lim", dest="dafbyaf_lim",
                  help="the prior distribution will be calculated from -dafbyaf_lim to -dafbyaf_lim", default=1)
  parser.add_option("--num-samples", dest="N_sampl",
                  help="number of prior samples to be used for computing the histogram", default=1e7)
  parser.add_option("--num-threads", dest="num_threads",
                  help="number of threads to be used", default=1)
  parser.add_option("--num-bins", dest="N_bins",
                  help="number of bins to be used for computing the histogram", default=401)
  (options, args) = parser.parse_args()

  fit_formula = options.fit_formula
  spin_angle_dist = options.spin_angle_dist
  uniform_in = options.uniform_in
  comp_mass_min = float(options.comp_mass_min)
  comp_mass_max = float(options.comp_mass_max)
  comp_spin_min = float(options.comp_spin_min)
  comp_spin_max = float(options.comp_spin_max)
  Mf_min = float(options.Mf_min)
  Mf_max = float(options.Mf_max)
  af_min = float(options.af_min)
  af_max = float(options.af_max)
  N_sampl = int(options.N_sampl)
  N_bins = int(options.N_bins)
  dMfbyMf_lim = float(options.dMfbyMf_lim)
  dafbyaf_lim = float(options.dafbyaf_lim)
  num_threads = int(options.num_threads)

  print '... N_sampl = %e N_bins = %d' %(N_sampl, N_bins)

  Mf_lim = max(abs(Mf_min), abs(Mf_max))
  af_lim = max(abs(af_min), abs(af_max))

  Mf_bins = np.linspace(-Mf_lim, Mf_lim, N_bins)
  af_bins = np.linspace(-af_lim, af_lim, N_bins)

  dMf = np.mean(np.diff(Mf_bins))
  daf = np.mean(np.diff(af_bins))

  Mf_intp = (Mf_bins[:-1] + Mf_bins[1:])/2.
  af_intp = (af_bins[:-1] + af_bins[1:])/2.

  if uniform_in == "Mfaf":

	print "... uniform in : final parameters"	

  	Mf_i, Mf_r = np.random.uniform(Mf_min, Mf_max, N_sampl), np.random.uniform(Mf_min, Mf_max, N_sampl)
  	af_i, af_r= np.random.uniform(af_min, af_max, N_sampl), np.random.uniform(af_min, af_max, N_sampl)

  elif uniform_in == "m1m2a1a2":

	print "... uniform in : initial parameters"	


	# generate random samples of the prior uniform in component masses 
  	m1_i = np.random.uniform(comp_mass_min, comp_mass_max, N_sampl)
  	m2_i = np.random.uniform(comp_mass_min, comp_mass_max, N_sampl)
  	m1_r = np.random.uniform(comp_mass_min, comp_mass_max, N_sampl)
  	m2_r = np.random.uniform(comp_mass_min, comp_mass_max, N_sampl)

	chi1_i = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)
	chi2_i = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)
	chi1_r = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)
        chi2_r = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)

	if spin_angle_dist == "isotropic":	
	  
	  print "... spin angle dist: isotropic"	

	  tilt1_i = np.arccos(np.random.uniform(-1, 1, N_sampl))
	  tilt2_i = np.arccos(np.random.uniform(-1, 1, N_sampl))
	  phi12_i = np.random.uniform(-1, 1, N_sampl)*np.pi

	  tilt1_r = np.arccos(np.random.uniform(-1, 1, N_sampl))
          tilt2_r = np.arccos(np.random.uniform(-1, 1, N_sampl))
          phi12_r = np.random.uniform(-1, 1, N_sampl)*np.pi

	  chi1z_i = chi1_i * np.cos(tilt1_i)
          chi2z_i = chi2_i * np.cos(tilt2_i)

          chi1z_r = chi1_r * np.cos(tilt1_r)
          chi2z_r = chi2_r * np.cos(tilt2_r)

	elif spin_angle_dist == "aligned":
  
	  print "... spin angle dist: aligned"	

          phi12_i = np.zeros(N_sampl)
          phi12_r = np.zeros(N_sampl)

	  chi1z_i = chi1_i * np.random.choice([-1,1],N_sampl)
	  chi2z_i = chi2_i * np.random.choice([-1,1],N_sampl)
          chi1z_r = chi1_r * np.random.choice([-1,1],N_sampl)
          chi2z_r = chi2_r * np.random.choice([-1,1],N_sampl)
	
	Mf_i, af_i = tgr.calc_final_mass_spin(m1_i, m2_i, chi1_i, chi2_i, chi1z_i, chi2z_i, phi12_i, fit_formula)
	Mf_r, af_r = tgr.calc_final_mass_spin(m1_r, m2_r, chi1_r, chi2_r, chi1z_r, chi2z_r, phi12_r, fit_formula)

  # compute the 2D posterior distributions for the inspiral, ringodwn and IMR analyses 
  P_Mfaf_i, Mf_bins, af_bins = np.histogram2d(Mf_i, af_i, bins=(Mf_bins, af_bins), normed=True)
  P_Mfaf_r, Mf_bins, af_bins = np.histogram2d(Mf_r, af_r, bins=(Mf_bins, af_bins), normed=True)

  # transpose to go from (X,Y) indexing returned by np.histogram2d() to array (i,j) indexing for further
  # computations. From now onwards, different rows (i) correspond to different values of Mf and different 
  # columns (j) correspond to different values of af 
  P_Mfaf_i = P_Mfaf_i.T
  P_Mfaf_r = P_Mfaf_r.T

  print '... created posteriors in P_Mfaf_i, P_Mfaf_r priors'

  # compute interpolation objects for the (delta_Mf, delta_af) posterior and (mean_Mf, mean_af) posterior 
  P_Mfaf_i_interp_object = interpolate.interp2d(Mf_intp, af_intp, P_Mfaf_i, fill_value=0., bounds_error=False)
  P_Mfaf_r_interp_object = interpolate.interp2d(Mf_intp, af_intp, P_Mfaf_r, fill_value=0., bounds_error=False)


  # defining limits of delta_Mf/Mf and delta_af/af.
  dMfbyMf_vec = np.linspace(-dMfbyMf_lim, dMfbyMf_lim, N_bins)
  dafbyaf_vec = np.linspace(-dafbyaf_lim, dafbyaf_lim, N_bins)

  # compute the P(dMf/Mf, daf/af) by evaluating the integral 
  diff_dMfbyMf = np.mean(np.diff(dMfbyMf_vec))
  diff_dafbyaf = np.mean(np.diff(dafbyaf_vec))
  P_dMfbyMfmean_dafbyafmean_pr = np.zeros(shape=(N_bins,N_bins))

  for i, v2 in enumerate(dafbyaf_vec):
    for j, v1 in enumerate(dMfbyMf_vec):
      P_dMfbyMfmean_dafbyafmean_pr[i,j] = tgr.calc_sum(Mf_intp, af_intp, v1, v2, P_Mfaf_i_interp_object, P_Mfaf_r_interp_object)

  # normalization
  P_dMfbyMfmean_dafbyafmean_pr /= np.sum(P_dMfbyMfmean_dafbyafmean_pr) * diff_dMfbyMf * diff_dafbyaf  
  
  print '... computed (delta_Mf/Mfmean, delta_af/afmean) prior'

  # create an interpolation object and save it 
  
  if uniform_in == "Mfaf":
  	outfile = 'Prior_normbymean_Mfmin_%.1f_Mfmax_%.1f_afmin%.1f_afmax%.1f_dMfbyMflim%.1f_dafbyaflim%.1f_%dbins_%s'%(Mf_min, Mf_max, af_min, af_max, dMfbyMf_lim, dafbyaf_lim, N_bins, uniform_in)
  elif uniform_in == "m1m2a1a2":
	outfile = 'Prior_normbymean_%s_m12min%.2f_m12max%.2f_a12min%.2f_a12max%.2f_epsilon%.1f_sigma%.1f_%dbins_%s'%(fit_formula, comp_mass_min, comp_mass_max,  comp_spin_min, comp_spin_max, dMfbyMf_lim, dafbyaf_lim, N_bins, uniform_in)

  P_dMfbyMfmean_dafbyafmean_pr_interp_obj = interpolate.interp2d(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMfmean_dafbyafmean_pr, fill_value=0., bounds_error=False)
  f = gzip.open(outfile+".pklz",'wb')
  pickle.dump(P_dMfbyMfmean_dafbyafmean_pr_interp_obj, f)

  print '... saved the interpolation object.'

  # read the interpolation object, reconstruct the data from the interpolation object 
  f = gzip.open(outfile+".pklz",'rb')
  P_dMfbyMfmean_dafbyafmean_pr_interp_obj = pickle.load(f)
  P_dMfbyMfmean_dafbyafmean_pr_interp = P_dMfbyMfmean_dafbyafmean_pr_interp_obj(dMfbyMf_vec, dafbyaf_vec)

  # difference between the original and interpolated data 
  interp_err = abs(P_dMfbyMfmean_dafbyafmean_pr - P_dMfbyMfmean_dafbyafmean_pr_interp)

  print '... maximum difference between the original and interpolated data is %e' %np.max(interp_err)

  # Random Sampling
  P_dMfbyMfmean_dafbyafmean_pr_rand, dMfbyMf_vec, dafbyaf_vec = np.histogram2d(2.*(Mf_i - Mf_r)/(Mf_i + Mf_r), 2*(af_i - af_r)/(af_i + af_r), bins=(dMfbyMf_vec, dafbyaf_vec), normed=True)

  P_dMfbyMfmean_dafbyafmean_pr_rand = P_dMfbyMfmean_dafbyafmean_pr_rand.T
  dMfbyMf_vec_intp = (dMfbyMf_vec[:-1] + dMfbyMf_vec[1:])/2.
  dafbyaf_vec_intp = (dafbyaf_vec[:-1] + dafbyaf_vec[1:])/2.


  plt.figure(figsize=(15,4))

  conf_v1v2 = confidence(P_dMfbyMfmean_dafbyafmean_pr_rand)
  s1_v1v2 = conf_v1v2.height_from_level(0.68)
  s2_v1v2 = conf_v1v2.height_from_level(0.95)

  plt.subplot(141)
  #plt.pcolormesh(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMfmean_dafbyafmean_pr_rand, cmap='YlOrBr')
  plt.scatter(2.*(Mf_i - Mf_r)/(Mf_i + Mf_r), 2*(af_i - af_r)/(af_i + af_r), marker='.', color='k', alpha=0.1)
  #plt.colorbar()
  #plt.clim(0, 1.)
  plt.contour(dMfbyMf_vec_intp, dafbyaf_vec_intp, P_dMfbyMfmean_dafbyafmean_pr_rand, levels=(s2_v1v2,s1_v1v2))
  plt.xlabel('$\Delta M_f / M_f$')
  plt.ylabel('$\Delta \chi _f / \chi _f$')
  plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
  plt.ylim(-dafbyaf_lim, dafbyaf_lim)
  plt.grid()
  plt.title('Random Sampling')


  conf_v1v2 = confidence(P_dMfbyMfmean_dafbyafmean_pr)
  s1_v1v2 = conf_v1v2.height_from_level(0.68)
  s2_v1v2 = conf_v1v2.height_from_level(0.95)

  plt.subplot(142)
  plt.pcolormesh(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMfmean_dafbyafmean_pr, cmap='YlOrBr')
  plt.colorbar()
  plt.clim(0, 1.)
  plt.contour(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMfmean_dafbyafmean_pr, levels=(s2_v1v2,s1_v1v2))
  plt.xlabel('$\Delta M_f / M_f$')
  plt.ylabel('$\Delta \chi _f / \chi _f$')
  plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
  plt.ylim(-dafbyaf_lim, dafbyaf_lim)
  plt.grid()
  plt.title('Original')

  conf_v1v2 = confidence(P_dMfbyMfmean_dafbyafmean_pr_interp)
  s1_v1v2 = conf_v1v2.height_from_level(0.68)
  s2_v1v2 = conf_v1v2.height_from_level(0.95)

  plt.subplot(143)
  plt.pcolormesh(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMfmean_dafbyafmean_pr_interp, cmap='YlOrBr')
  plt.colorbar()
  plt.clim(0, 1.)
  plt.contour(dMfbyMf_vec, dafbyaf_vec, P_dMfbyMfmean_dafbyafmean_pr_interp, levels=(s2_v1v2,s1_v1v2))
  plt.xlabel('$\Delta M_f / M_f$')
  plt.ylabel('$\Delta \chi _f / \chi _f$')
  plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
  plt.ylim(-dafbyaf_lim, dafbyaf_lim)
  plt.grid()
  plt.title('Interpolation')

  plt.subplot(144)
  plt.pcolormesh(dMfbyMf_vec, dafbyaf_vec, interp_err, cmap='YlOrBr')
  plt.xlabel('$\Delta M_f / M_f$')
  plt.ylabel('$\Delta \chi _f / \chi _f$')
  plt.xlim(-dMfbyMf_lim, dMfbyMf_lim)
  plt.ylim(-dafbyaf_lim, dafbyaf_lim)
  plt.grid()
  plt.clim(np.min(interp_err)-1e-16,np.max(interp_err)+1e-16)
  plt.colorbar()
  plt.title('Difference')

  plt.savefig(outfile+'.png', dpi=300)
  print '... saved fig (time taken: %f secs)' %(time.time()-start_time)
