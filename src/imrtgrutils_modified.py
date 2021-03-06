"""
Core python module for the IMR consistency test 

P. Ajith, Abhirup Ghosh, 2015-09-15

$Id:$
"""

import numpy as np 
import scipy.ndimage.filters as filter
from multiprocessing import Pool
from functools import partial
import time 
import nrutils as nr

""" calculate the mass and spin of the final black hole using an NR-inspired fitting formula """
def calc_final_mass_spin(m1, m2, chi1, chi2, fit_formula):
  
  if fit_formula == 'nospin_Pan2011':
    af = nr.bbh_final_spin_non_spinning(m1, m2)
    mf = nr.bbh_final_mass_non_spinning(m1, m2)
  elif fit_formula == 'nonprecspin_Healy2014':
    af = nr.bbh_final_spin_non_precessing(m1, m2, chi1, chi2)
    mf = nr.bbh_final_mass_non_precessing(m1, m2, chi1, chi2, af)
  else:
    print '# unknown spin fit formula'
    exit()
  return mf, af

""" compute the integrant of P(dMf/Mf, daf/af). """
def P_integrant(af, Mf, v1, v2, P_dMfdaf_interp_object, P_Mfaf_imr_interp_object):

  Mf_mat, af_mat = np.meshgrid(Mf, af)

  # Create dMf and daf vectors corresponding to the given v1 and v2. These vectors have to be 
  # monotonically increasing in order to evaluate the interpolated prob densities. Hence, for 
  # v1, v2 < 0, flip them, evaluate the prob density (in column or row) and flip it back 
  dMf = v1*Mf
  daf = v2*np.ones(len(af))

  if v1 < 0.:
    dMf = np.flipud(dMf)
  if v2 < 0.:
    daf = np.flipud(daf)
  P_delta = P_dMfdaf_interp_object(dMf, daf)

  if v1 < 0.:
    P_delta = np.fliplr(P_delta)
  if v2 < 0.:
    P_delta = np.flipud(P_delta)

  P_imr = P_Mfaf_imr_interp_object(Mf, af)
  return P_imr*P_delta*abs(Mf_mat), P_imr, P_delta

""" compute P(dMf/Mf, daf/af). """
def calc_sum(Mf, af, v1, v2, P_dMfdaf_interp_object, P_Mfaf_imr_interp_object):

  Pintg, P_imr, P_delta = P_integrant(af, Mf, v1, v2, P_dMfdaf_interp_object, P_Mfaf_imr_interp_object)
  return np.sum(Pintg)

""" gaussian filter of histogram """
def gf(P):
  return filter.gaussian_filter(P, sigma=2.0)

""" generate prior samples in (Mf, af) assuming uniform prior in component masses """
def calc_Mfaf_prior_samples(comp_mass_prior_min, comp_mass_prior_max, comp_spin_min, comp_spin_max, Mf_bins, af_bins, fit_formula, N_sampl, thread):

  # generate random samples of the prior uniform in component masses 
  m1_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, N_sampl)
  m2_pr = np.random.uniform(comp_mass_prior_min, comp_mass_prior_max, N_sampl)
  chi1_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)
  chi2_pr = np.random.uniform(comp_spin_min, comp_spin_max, N_sampl)

  # return the corrresponding samples of prior in Mf, af 
  return calc_final_mass_spin(m1_pr, m2_pr, chi1_pr, chi2_pr, fit_formula)

""" calculate the prior distribution in (Mf, af) assuming uniform prior in component masses """
def calc_Mfaf_prior(comp_mass_prior_min, comp_mass_prior_max, comp_spin_min, comp_spin_max, Mf_bins, af_bins, fit_formula, N_sampl, num_threads):

  # generate samples in Mf and af corresponding to uniform samples in m1, m2. Parallelise the calculation over multiple threads 
  thread_vec = range(num_threads)
  N_sampl = int(N_sampl/num_threads)
  p = Pool(num_threads)
  func = partial(calc_Mfaf_prior_samples, comp_mass_prior_min, comp_mass_prior_max, comp_spin_min, comp_spin_max, Mf_bins, af_bins, fit_formula, N_sampl)
  final_params = p.map(func, thread_vec)
  p.close()
  p.join()
  
  final_params = np.array(final_params)
  final_params.reshape(N_sampl*num_threads, 2)
  Mf_pr = np.ravel(final_params[:,0])
  af_pr = np.ravel(final_params[:,1])

  # compute the 2D prior distribution in Mf and af 
  P_Mfaf_pr, Mf_bins, af_bins = np.histogram2d(Mf_pr, af_pr, bins=(Mf_bins, af_bins), normed=True)

  return P_Mfaf_pr.T
